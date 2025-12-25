# Aero Standard Paths - Final Implementation

## Coordinate Conventions
- The arena is in NED coordinate system in meters
- Entry point is always heading south at x=0, y=0, z=-25
- +x is north (12 o'clock when viewed from above)
- +y is east (3 o'clock)
- +z is down (increasing depth)
- Right turn when heading south = turning to the west (clockwise from above)
- Pitch up when heading south = going toward more negative z (upward)

## Path Details
Path details are in terms of destination points or maneuvers from the current position:
- A right turn is a turn in the constant z plane that moves the heading clockwise
- A pitch up is a loop that moves craft nose up (toward more -z)
- The arena is generally 100m in each direction but we stay inside as described in other path generators

## Final Path Set (6 Paths - IMPLEMENTED)

### Path 0: Straight and Level ✓
Racetrack pattern at constant altitude z=-25
- Head south for 20 meters dest=[-20,0,-25]
- Execute a right 180° turn with 20m radius dest=[-20,-40,-25]
- Head north for 40 meters dest=[20,-40,-25]
- Execute a right 180° turn with 20m radius dest=[20,0,-25]
- Head back to origin for 20 meters dest=[0,0,-25]
- **Key**: All segments maintain z=-25 (no altitude change)

### Path 1: Spiral Climb ✓
1.5 circles climbing maneuver
- Head straight south for 20 meters dest=[-20,0,-25]
- Execute a 540° turn (1.5 circles) clockwise from above, with 20m radius, smoothly climbing to z=-75
- Head north for 40 meters at z=-75
- **Key**: Clockwise from above (right turn from pilot perspective)

### Path 2: Horizontal Figure Eight ✓
Two horizontal loops (left then right)
- Lead-in: 1 second straight south
- Left horizontal loop (r=20m)
- Right horizontal loop (r=20m)
- **Key**: Reuses longSequential logic

### Path 3: 45 Degree Angled Loop ✓
Diagonal loop tilted 45° combining horizontal and vertical motion
- Loop rotated -45° around x-axis
- Entry at bottom of loop at origin (0, 0, -25)
- Loop goes UP (into negative Z, higher altitude)
- Continues southward from entry
- **Key**: Same as FortyFiveLoop in computedPaths but with bottom entry point at origin

### Path 4: High Perch with Split-S Reversal
Complex multi-segment maneuver combining climb, perch, and reversal
1. **Initial South**: Head south for 20m to [-20, 0, -25]
2. **Climbing Left Turn**: 180° left turn (clockwise from top view) with 20m radius while climbing 20m to z=-45
   - Left turn from pilot heading south = clockwise when viewed from -z looking +z
   - Ends at approximately [-20, -40, -45]
3. **Diagonal Climb**: Straight climb heading north, ascending from z=-45 to z=-65 over ~40m horizontal distance
   - Ends at approximately [20, -40, -65] (high perch position)
4. **Hard Right Turn**: Tight 150° right turn with 5m radius
   - From north heading, turns 150° clockwise to southwest heading
5. **Straight Segment**: Fly 30m in new heading (150° from north = southwest)
6. **Continue to x=-20**: Continue in same heading until reaching x=-20
7. **Pitch Down at x=-20**: Execute 180° pitch-down loop (course reversal in +z direction)
   - Split-S maneuver: half loop downward, reversing course
   - Inverted entry, rolling upright while diving
8. **Final North Run**: Head north to x=40
- **Key**: Combines climbing spiral, perch maneuver, and split-S reversal

### Path 5: Seeded Random B ✓
Reproducible random path
- Generated using random path generator with FIXED seed=67890
- Provides variety while maintaining reproducibility
- Seed configurable via `RandomPathSeedB` in autoc.ini

## Removed Paths (Too Complex)
- ~~Cuban Eight~~ - Complex pitch loops, removed for simplicity
- ~~Split S / Immelmann~~ - Complex maneuver, removed for simplicity
- ~~Rapid Course Reversal~~ - Instant turns not smooth, removed

## Implementation Details

### Technique
- Entry: Virtually heading south at x=0, y=0, z=-25
- Course stays within arena limits
- New course type: `aeroStandard` alongside `longSequential`
- Base speed: Same as other path types
- Sample rate: Roughly the same as other generators

### Key Fixes Applied
1. **Z Discontinuity Fix**: Fixed bug where horizontal turns were doubling Z coordinate
2. **Path 3 Direction Fix**: Loop now heads southward (negative x) from entry
3. **Path 3 Rotation Fix**: Loop goes UP (negative z) instead of down into ground
4. **Path 3 Position Fix**: Bottom of loop placed at origin (0, 0, -25)

### Configuration
```ini
SimNumPathsPerGeneration = 6
PathGeneratorMethod = aeroStandard
RandomPathSeedB = 67890
```

## Goals
- More complex set of training paths
- Much more all-attitude maneuvering
- Later: Add entry variations and craft variations
- Target: Closer to zero-shot generalization in variations project

## Code Structure
- **Enum**: `AeroStandardPathType` with 6 path types
- **Class**: `GenerateAeroStandard` in pathgen.h
- **Helper Methods**:
  - `addStraightSegment()` - Linear segments
  - `addHorizontalTurn()` - Constant-Z turns
  - `addSpiralTurn()` - Turns with altitude change
  - `addHorizontalLoop()` - 360° circles

## Continuity
- All paths use rabbit/turtle mode (pen never lifts)
- Each segment starts where previous segment ended
- No discontinuous jumps

## Training Use
The system generates 6 diverse, reproducible paths per generation for GP evolution.
