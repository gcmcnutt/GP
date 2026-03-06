# FLASH_SPEEDUP: WebBLE Download Performance

## Status: IMPLEMENTED - Testing Simple Request-Response

## Problem

Current WebBLE download speed: **1-2 KB/sec**
Typical flight log: **200-500 KB** (3-5 min at 10Hz, ~163 bytes/message)
Download time: **2-8 minutes** (unacceptable for field use)

### Baseline Measurement
- **17.9 KB in 18 seconds = 0.99 KB/s**

## Key Insight: Walk or Chew Gum

The system enforces **mutually exclusive modes**:
- **Armed (logging)**: Flash writes active, BLE disabled, downloads blocked
- **Disarmed (idle)**: No writes, BLE enabled, downloads allowed

Code enforcement:
- `flashLoggerStartDownload()` rejects if `!loggingSuspended` (flash_logger.cpp:481-484)
- `blueToothSetEnabled(false)` during armed state

**This means read/write conflicts are impossible.** Much of the current defensive code was added during debugging when the actual root cause was simultaneous read/write attempts. Now that the architecture enforces mutual exclusion, these safeguards are unnecessary overhead.

## Root Cause Analysis

### Hardware Specs

**Flash chip: P25Q16H**
- Max speed: 55 MHz
- Current config: 4 MHz (`NRF_QSPI_FREQ_32MDIV8`)

**nRF52840 QSPI options:**
| Divider | Frequency |
|---------|-----------|
| 32MDIV1 | 32 MHz |
| 32MDIV2 | 16 MHz |
| 32MDIV4 | 8 MHz |
| 32MDIV8 | 4 MHz (current) |

### Bottlenecks (All Removable Due to Walk/Chew Gum)

1. **Pre-download file scan** - `computeFilteredSize()` reads ENTIRE file before transfer
   - *Unnecessary*: When disarmed, file size is stable - just use raw endAddr-startAddr
2. **Per-byte null filtering** - Loops through every byte during read
   - *Unnecessary*: Nulls exist only as word-alignment padding at buffer boundaries
   - Alternative: Track actual payload bytes written, not padded size
3. **Sync calls during LIST** - `flashLoggerGetFileCount/Size()` each call 250ms sync
   - *Unnecessary*: When disarmed, nothing is being written - no flush needed
4. **10ms delays in LIST** - Hardcoded `delay(10)` per file notification
   - *Overcautious*: Defensive timing from debug phase
5. **Request-response per chunk** - Each 128 bytes requires full BLE round-trip (~20-30ms)
   - *Inefficient*: Streaming mode would eliminate round-trips
6. **Small chunk size** - 128 bytes vs BLE MTU of 240-512 bytes
   - *Conservative*: Comment says "smaller for reliability" - legacy from debug phase
7. **QSPI at 4 MHz** - Flash supports 55 MHz, nRF52840 supports 32 MHz
   - *Conservative*: Safe to increase since no concurrent operations during download

## Simplified Implementation Plan

Given the walk/chew-gum architecture, we can aggressively remove defensive code.

### Phase 1: Remove Defensive Overhead (Firmware Only)

**No protocol changes, just cleanup:**

```cpp
// flash_logger.cpp L145 - Increase QSPI frequency
.sck_freq = NRF_QSPI_FREQ_32MDIV2  // 16 MHz (safe) or MDIV1 for 32 MHz

// flash_logger.h - Larger chunks (matches BLE MTU)
#define FLASH_LOGGER_CHUNK_SIZE 240

// bluetooth.cpp L108, L117 - Reduce LIST delays
delay(2);  // Was 10ms

// flash_logger.cpp flashLoggerStartDownload() - Skip pre-scan
// REMOVE: uint32_t filteredLength = computeFilteredSize(...)
// INSTEAD: Just use raw size
downloadFilteredSize = endAddr - startAddr;

// flash_logger.cpp flashLoggerReadChunk() - Remove null filtering
// REMOVE: The byte-by-byte null check loop
// INSTEAD: Just copy raw bytes (nulls are only alignment padding at boundaries)

// flash_logger.cpp flashLoggerGetFileCount/Size() - Remove sync calls
// REMOVE: flashLoggerSyncBlocking(250ms) at start
// When disarmed, nothing is being written - sync is pointless
```

**Expected improvement: 1-2 KB/s → 15-20 KB/s**

### Phase 2: Windowed Flow Control (ABANDONED)

**Attempted but abandoned due to reliability issues.**

Pure streaming overflows BLE stack buffers. Tried windowed flow control:
- Firmware sends N chunks (window), then sends "WINDOW" status
- Browser sends "ACK" command when ready
- Firmware sends next window

**Problems encountered:**
- Even with WINDOW_SIZE=4 and delays (3ms inter-chunk, 20ms before WINDOW), transfers would intermittently wedge
- Out-of-order notifications between data and status characteristics
- CRC mismatches indicating lost data
- The ArduinoBLE stack doesn't provide reliable flow control primitives

**Research finding:** "Many stacks drop data within the software stack due to memory constraints and queue limitations" — [Memfault BLE Throughput](https://interrupt.memfault.com/blog/ble-throughput-primer)

**Conclusion:** Simple request-response is the most reliable approach. The ~1.64 KB/s achieved in Phase 1 is acceptable for typical flight logs (17.9 KB in ~11 sec). Focus should be on reducing the per-request overhead rather than batching.

### Current Approach: Simple Request-Response (Zero Delays)

Since we have pure sync (browser requests chunk → firmware sends chunk → browser requests next), no delays are needed:

```cpp
// bluetooth.cpp - Simple one-chunk-per-request
static void transferNextChunk() {
  int bytesRead = flashLoggerReadChunk(transferBuffer, FLASH_LOGGER_CHUNK_SIZE);
  if (bytesRead < 0) { /* error */ return; }
  if (bytesRead == 0) { /* send DOWNLOAD_DONE with CRC */ return; }

  downloadCrc = crc32Update(downloadCrc, transferBuffer, bytesRead);
  downloadTotalBytes += bytesRead;
  dataCharacteristic.writeValue(transferBuffer, bytesRead);
  // No delay - wait for browser NEXT
}
```

Browser requests each chunk:
```javascript
function handleDataChunk(event) {
  // accumulate data
  requestNextChunk();  // immediately request next
}

async function requestNextChunk() {
  await controlChar.writeValue(textEncoder.encode('NEXT'));
}
```

**Benefits of simple request-response:**
- Reliable: pure sync, no race conditions
- Simple: no windowing state to manage
- Debuggable: clear request-response flow

### Phase 3: Compression (Future, If Needed)

Skip for now. Raw transfer at 30-50 KB/s means:
- 300KB file → 6-10 seconds (acceptable for field use)

If compression is ever needed:
- **Heatshrink**: 50 bytes RAM, Arduino library available
- **Browser**: Native DecompressionStream API (Chrome 80+, zero overhead)

### Why This Simplification Works

The original defensive code was written when the system had potential for concurrent read/write. Now that architecture enforces mutual exclusion:

| Concern | Original Defense | Why It's Safe to Remove |
|---------|------------------|------------------------|
| File size changing during download | Pre-scan with `computeFilteredSize()` | Disarmed = no writes = stable size |
| Incomplete flush | `flashLoggerSyncBlocking()` in LIST | Disarmed = nothing to flush |
| Data corruption | Per-byte null filtering | No concurrent writes to corrupt data |
| BLE timing issues | 10ms delays, small chunks | No flash contention during download |

## Files to Modify

### Phase 1 (Cleanup)
- `src/flash_logger.cpp`:
  - Line 145: `NRF_QSPI_FREQ_32MDIV8` → `NRF_QSPI_FREQ_32MDIV2` (16 MHz)
  - Line 532-544: Remove `computeFilteredSize()` call, use raw size
  - Lines 586-592: Remove null-filtering loop in `flashLoggerReadChunk()`
  - Lines 433, 463: Remove `flashLoggerSyncBlocking()` calls in file info functions
- `include/flash_logger.h`:
  - Line 7: `FLASH_LOGGER_CHUNK_SIZE 128` → `240`
- `src/bluetooth.cpp`:
  - Lines 108, 117: `delay(10)` → `delay(2)`

### Phase 1 (Browser Cleanup)
- `web/flight_logger.html`:
  - Remove trailing zero trim loop (lines 366-372) - firmware handles null filtering
  - Throttle `updateProgress()` - call every 10 chunks or 100ms, not every chunk
  - Cache `TextEncoder` instance instead of recreating per command
  - Add CRC32 verification (see Checksum section)

### Phase 2 (Windowed) - ABANDONED
See "Phase 2: Windowed Flow Control (ABANDONED)" section above for details on why this approach was not viable.

## Checksum Implementation (Required)

Full-file CRC32 is simple and sufficient given walk/chew-gum architecture.

### Firmware Side (flash_logger.cpp / bluetooth.cpp)

```cpp
// During download, accumulate CRC32
#include <CRC32.h>  // Arduino CRC32 library

static CRC32 downloadCrc;

// In flashLoggerStartDownload():
downloadCrc.reset();

// In flashLoggerReadChunk(), after reading each chunk:
downloadCrc.update(buffer, bytesRead);

// In transferNextChunk() when bytesRead == 0 (complete):
char statusMsg[64];
snprintf(statusMsg, sizeof(statusMsg), "DOWNLOAD_DONE:%08lX:%lu",
         downloadCrc.finalize(), totalBytesSent);
statusCharacteristic.writeValue(statusMsg);
```

### Browser Side (flight_logger.html)

```javascript
// Add CRC32 function (fast, ~20 lines, no library needed)
function crc32(data) {
  let crc = 0xFFFFFFFF;
  for (let i = 0; i < data.length; i++) {
    crc ^= data[i];
    for (let j = 0; j < 8; j++) {
      crc = (crc >>> 1) ^ (crc & 1 ? 0xEDB88320 : 0);
    }
  }
  return (crc ^ 0xFFFFFFFF) >>> 0;
}

// In handleStatusChange(), parse DOWNLOAD_DONE:
} else if (value.startsWith('DOWNLOAD_DONE:')) {
  const parts = value.split(':');
  const expectedCrc = parseInt(parts[1], 16);
  const expectedSize = parseInt(parts[2]);

  const actualCrc = crc32(downloadBuffer.subarray(0, downloadOffset));

  if (actualCrc !== expectedCrc) {
    log(`ERROR: CRC mismatch! Expected ${expectedCrc.toString(16)}, got ${actualCrc.toString(16)}`);
    hideProgress();
    return;
  }

  if (downloadOffset !== expectedSize) {
    log(`WARNING: Size mismatch. Expected ${expectedSize}, got ${downloadOffset}`);
  }

  saveDownload();
  hideProgress();
}
```

### Why Full-File (Not Per-Block)

| Approach | Pro | Con |
|----------|-----|-----|
| Per-block CRC | Detect corruption early, retry | Protocol complexity, overhead |
| Full-file CRC | Simple, no per-chunk overhead | Only detects at end |

With walk/chew-gum, corruption is unlikely (no concurrent write). Full-file CRC catches flash read errors or BLE transmission errors with minimal complexity.

## Performance Results

| Metric | Baseline | Phase 1 (Final) | Improvement |
|--------|----------|-----------------|-------------|
| Speed | **0.99 KB/s** | **1.77 KB/s** | 1.79x |
| 17.9KB test file | 18 sec | ~10 sec | 1.8x |
| 171.4KB test file | ~173 sec | **97 sec** | 1.78x |
| 300KB flight log | ~5 min | ~2.8 min | 1.8x |
| Protocol | Request/response | Request/response | Zero delays |

**Key findings:**
1. Request-response overhead dominates (~135ms per round-trip at 240 bytes/chunk)
2. Windowed flow control attempted but unreliable (wedging, data loss)
3. Simple request-response with Phase 1 cleanup (QSPI 16MHz, 240-byte chunks, no pre-scan) provides best reliability/speed tradeoff
4. CRC32 verification ensures data integrity

## Phase 3: PRN-Style Pipelining (Nordic DFU Inspired)

Based on Nordic's proven DFU protocol which uses **Packet Receipt Notification (PRN)** for reliable high-speed transfers.

### Concept: Windowed Requests with Periodic CRC Checkpoints

Instead of firmware-driven windowing (which caused wedging), the **browser drives pipelining**:

```
Browser                              Firmware
   |---NEXT (fire-and-forget)-------->|
   |---NEXT (don't wait)------------->|
   |---NEXT (pipeline ahead)--------->|
   |<---------[chunk 0]---------------|
   |---NEXT (keep pipeline full)----->|
   |<---------[chunk 1]---------------|
   |---NEXT---------------------------|
   |<---------[chunk 2]---------------|
   ... (after PRN_INTERVAL chunks) ...
   |---CRC (checkpoint request)------>|
   |<----CRC:offset:crc---------------|  (validate, continue or retry)
   |---NEXT-------------------------->|
   ...
```

**Key insight**: Browser controls pacing. If BLE stack drops data, CRC checkpoint catches it.

### Protocol

1. **NEXT** - Request next chunk (use `writeValueWithoutResponse` for speed)
2. **CRC** - Request current offset + cumulative CRC (use `writeValue` to ensure delivery)
3. Firmware responds to CRC with `CRC:<offset>:<crc32>` on status characteristic

### Browser Implementation

```javascript
const PRN_INTERVAL = 8;      // CRC checkpoint every N chunks
const PIPELINE_DEPTH = 4;    // Outstanding NEXT requests
let outstandingRequests = 0;
let chunksSinceCheckpoint = 0;
let localCrc = 0xFFFFFFFF;
let lastValidatedOffset = 0;

function startPipeline() {
  outstandingRequests = 0;
  chunksSinceCheckpoint = 0;
  // Prime the pipeline
  for (let i = 0; i < PIPELINE_DEPTH && downloadOffset < expectedFileSize; i++) {
    sendNextRequest();
  }
}

function sendNextRequest() {
  if (downloadOffset + (outstandingRequests * CHUNK_SIZE) < expectedFileSize) {
    controlChar.writeValueWithoutResponse(textEncoder.encode('NEXT'));
    outstandingRequests++;
  }
}

function handleDataChunk(event) {
  // Accumulate data
  const chunk = new Uint8Array(event.target.value.buffer);
  downloadBuffer.set(chunk, downloadOffset);
  downloadOffset += chunk.length;
  outstandingRequests--;

  // Update local CRC
  localCrc = crc32Update(localCrc, chunk);
  chunksSinceCheckpoint++;

  // Periodic checkpoint
  if (chunksSinceCheckpoint >= PRN_INTERVAL) {
    requestCrcCheckpoint();
  } else {
    sendNextRequest();  // Keep pipeline full
  }
}

async function requestCrcCheckpoint() {
  chunksSinceCheckpoint = 0;
  await controlChar.writeValue(textEncoder.encode('CRC'));
  // Wait for CRC response via status notification
}

function handleCrcResponse(offset, firmwareCrc) {
  const expectedCrc = (localCrc ^ 0xFFFFFFFF) >>> 0;
  if (offset !== downloadOffset || firmwareCrc !== expectedCrc) {
    console.error(`CRC mismatch at offset ${offset}: expected ${expectedCrc.toString(16)}, got ${firmwareCrc.toString(16)}`);
    // Could implement retry from lastValidatedOffset here
    return false;
  }
  lastValidatedOffset = offset;
  // Continue pipeline
  for (let i = 0; i < PIPELINE_DEPTH; i++) {
    sendNextRequest();
  }
  return true;
}
```

### Firmware Implementation

```cpp
// In processLoggerCommand():
} else if (strncmp(cmd, "CRC", 3) == 0) {
  // Respond with current download state for validation
  uint32_t currentCrc = downloadCrc ^ 0xFFFFFFFF;  // Finalize for comparison
  char msg[48];
  snprintf(msg, sizeof(msg), "CRC:%lu:%08lX",
           (unsigned long)downloadTotalBytes, (unsigned long)currentCrc);
  statusCharacteristic.writeValue(msg);
}
```

### Expected Performance

| Metric | Current | PRN Pipeline | Improvement |
|--------|---------|--------------|-------------|
| Requests/chunk | 1 | 0.125 (1 per 8) | 8x less overhead |
| Round-trip wait | Every chunk | Every 8 chunks | Hide latency |
| Speed | 1.77 KB/s | 4-7 KB/s | 2-4x |
| 171KB file | 97 sec | 25-40 sec | 2.5-4x |

### Experimental Results (Feb 2026)

**Prototype implemented but ABANDONED due to complexity and reliability issues.**

Actual results:
- 17.9KB file: 7.1 sec = **2.52 KB/s** (1.42x improvement, not 2-4x)
- 171KB file: **wedges frequently**, unreliable

Issues encountered:

1. **ArduinoBLE `written()` doesn't queue** - multiple rapid `writeValueWithoutResponse` calls overwrite each other. Only the last write is seen by firmware. Had to fall back to `writeValue` (with response), negating most pipelining benefit.

2. **State tracking complexity** - tracking `outstandingRequests`, `pendingRequestCount`, `awaitingDrain`, `awaitingCrcResponse`, `pendingCrcRequest`, `isSendingRequests` leads to race conditions. Chunks arriving during async `writeValue` calls cause count mismatches.

3. **Offset drift** - even with careful ordering, firmware offset can drift ahead of browser offset when chunks arrive during CRC request, causing checkpoint failures.

4. **Bursty but unreliable** - short files complete ~40% faster, but longer files wedge due to accumulated state machine bugs.

**Conclusion:** The 1.42x improvement doesn't justify the complexity. Simple request-response at 1.77 KB/s is reliable and acceptable. For future speedup, consider:
- Firmware-side command queuing (not possible with ArduinoBLE)
- Larger chunk sizes (limited by BLE MTU)
- Compression (Heatshrink on firmware, DecompressionStream on browser)

### References
- [Nordic DFU PRN](https://devzone.nordicsemi.com/f/nordic-q-a/45411/dfu-control-point---set-packet-receipt-notification)
- [Nordic SDK BLE Transport](https://infocenter.nordicsemi.com/topic/sdk_nrf5_v17.1.0/lib_dfu_transport_ble.html)
- [web-bluetooth-dfu](https://github.com/thegecko/web-bluetooth-dfu) - Working WebBLE implementation
- [writeValueWithoutResponse](https://developer.mozilla.org/en-US/docs/Web/API/BluetoothRemoteGATTCharacteristic/writeValueWithoutResponse)

## Browser-Side Cleanup Details

### Remove Trailing Zero Trim (lines 366-372)

```javascript
// REMOVE this loop in saveDownload():
let trimCount = 0;
while (downloadOffset > 0 &&
       downloadBuffer[downloadOffset - 1] === 0 &&
       trimCount < 512) {
  downloadOffset--;
  trimCount++;
}
```

If firmware keeps null filtering, this is redundant. If firmware removes null filtering, the CRC would include the nulls anyway, so trimming would cause CRC mismatch. Either way, remove it.

### Throttle Progress Updates

```javascript
// In handleDataChunk(), throttle DOM updates:
let lastProgressUpdate = 0;
const PROGRESS_THROTTLE_MS = 100;

function handleDataChunk(event) {
  // ... accumulate bytes ...

  const now = Date.now();
  if (now - lastProgressUpdate > PROGRESS_THROTTLE_MS) {
    updateProgress(...);
    lastProgressUpdate = now;
  }

  // ... request next chunk (Phase 1) or just return (Phase 2) ...
}
```

At 20 KB/s with 240-byte chunks, that's ~83 chunks/sec. DOM updates at 10 Hz instead of 83 Hz.

### Cache TextEncoder

```javascript
// At top of script:
const textEncoder = new TextEncoder();

// Then use textEncoder.encode() instead of new TextEncoder().encode()
```

Minor optimization but cleaner.

## Null Byte Handling (Important Nuance)

The current code filters null bytes for two reasons:
1. **Word-alignment padding**: Payloads are padded to 4-byte boundaries with nulls
2. **Unused flash regions**: Erased flash reads as 0xFF, but if a flight ends mid-buffer, remaining buffer bytes are 0x00

**Option A: Keep null filtering (simple)**
- Retain the null-filter loop but it's now faster because no concurrent read/write pressure
- Browser receives clean text without embedded nulls

**Option B: Track actual payload length (more work)**
- Store actual payload bytes written per flight (not padded size)
- Requires metadata schema change
- Avoids byte-by-byte loop entirely

**Recommendation**: Keep null filtering for Phase 1 (minimal risk). The overhead is small compared to pre-scan removal. Revisit if benchmarks show it's a bottleneck.

## Risk Mitigation

| Change | Risk | Mitigation |
|--------|------|------------|
| QSPI 4→16 MHz | Flash timing issues | Test with 8 MHz first (MDIV4) |
| Remove pre-scan | Browser expects wrong size | Size is now raw; browser handles it |
| Remove sync in LIST | Stale data reported | Only affects in-progress flight (disarmed = no in-progress) |
| Remove trailing zero trim | File has embedded nulls | CRC verification catches any mismatch |
| CRC mismatch | Bad download silently used | Browser refuses to save, shows error |

## References

- P25Q16H datasheet: https://files.seeedstudio.com/wiki/github_weiruanexample/Flash_P25Q16H-UXH-IR_Datasheet.pdf
- nRF52840 QSPI: https://docs.nordicsemi.com/bundle/ps_nrf52840/page/qspi.html
- Seeed QSPI usage: https://wiki.seeedstudio.com/xiao-ble-qspi-flash-usage/
- Heatshrink: https://github.com/atomicobject/heatshrink
- fflate: https://github.com/101arrowz/fflate
- Compression Streams API: https://developer.chrome.com/blog/compression-streams-api/
