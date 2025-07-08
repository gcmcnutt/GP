# Blackbox Integration for GP Renderer

## Overview

The GP renderer has been extended to support rendering real flight path data from iNav blackbox logs alongside the simulated flight paths. This allows visualization of actual flight data overlaid on the expected paths and simulated paths.

## Features

- **Conditional Blackbox Rendering**: Green flight path only appears when blackbox file is specified
- **Blackbox Data Extraction**: Parses iNav blackbox log files to extract position and attitude data
- **GPS to Local Conversion**: Converts GPS coordinates to local coordinate system centered on the initial position
- **Visual Overlay**: Renders blackbox data as a green tape/ribbon overlay on the simulation arena
- **Command Line Integration**: Simple command line option to load blackbox files
- **No Dependencies**: Renderer works normally without blackbox data, no errors or warnings

## Usage

### Basic Usage

```bash
# Run renderer with simulation data only
./renderer [simulation_key]

# Run renderer with simulation data and blackbox overlay
./renderer [simulation_key] [blackbox_file.TXT]
```

### Example

```bash
# Load specific simulation run with blackbox overlay
./renderer autoc-20250708-123456/ /path/to/blackbox_log_2025-07-08_143721.TXT

# Load latest simulation run with blackbox overlay
./renderer "" /path/to/blackbox_log_2025-07-08_143721.TXT

# Get help
./renderer --help
```

## Blackbox File Format

The extractor expects iNav blackbox log files in CSV format with the following structure:

- Header line starting with `loopIteration` 
- Data lines with comma-separated values including:
  - Timestamp (milliseconds)
  - GPS coordinates (latitude, longitude, altitude)
  - Attitude data (roll, pitch, yaw)
  - Control inputs (throttle, roll/pitch commands)

## Implementation Details

### Key Components

1. **BlackboxExtractor**: Core class for parsing blackbox files
2. **SimpleAircraftState**: Simplified aircraft state structure for blackbox data
3. **GPS Conversion**: Converts GPS coordinates to local NED coordinates
4. **Renderer Integration**: Adds blackbox data as a fourth visual layer

### Coordinate System

- **Input**: GPS coordinates (WGS84) + attitude in degrees
- **Output**: Local coordinates centered on first GPS position
- **Visualization**: X=East, Y=North, Z=Up (renderer coordinate system)

### Visual Representation

- **Expected Paths**: Red tubes (simulation target paths)
- **Actual Simulation**: Yellow/Orange tape (simulated aircraft flight)
- **Blackbox Data**: Green tape (real flight data)
- **Segments**: Blue lines connecting simulation states to expected paths

## Testing

A test program is included to verify blackbox extraction:

```bash
./test_blackbox_simple
```

## Future Enhancements

- Support for different blackbox field configurations
- Time synchronization between simulation and blackbox data
- Performance metrics comparison (simulation vs. actual)
- Multiple blackbox file support for flight comparisons
- Interactive controls for toggling data layers

## Dependencies

- Eigen3 (for vector/quaternion math)
- VTK (for visualization)
- Standard C++ libraries (fstream, regex, etc.)