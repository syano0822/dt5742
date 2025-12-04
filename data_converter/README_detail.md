# Waveform Processing Pipeline

Multi-stage waveform data converter and analyzer for CAEN DT5742 digitizer data with modular C++ architecture.

## Overview

This pipeline processes raw digitizer data through three stages:

1. **waveform_converter** (`src/convert_to_root.cpp`): Binary/ASCII → ROOT (raw waveforms with pedestal subtraction)
2. **waveform_analyzer** (`src/analyze_waveforms.cpp`): ROOT → ROOT (timing and amplitude feature extraction)
3. **hdf5_exporter** (`src/export_to_hdf5.cpp`): ROOT → HDF5 (export for Python/NumPy analysis)

## Quick Start

```bash
# 1. Setup environment
source /opt/root6/root_install/bin/thisroot.sh

# 2. Build
make

# 3. Configure
# Edit converter_config.json (unified configuration for all stages)

# 4. Run full pipeline
./run_full_pipeline.sh

# Or run stages individually (all use converter_config.json by default)
./convert_to_root
./analyze_waveforms
./export_to_hdf5 --mode raw
./export_to_hdf5 --mode analysis
```

## Project Structure

```
data_converter/
├── include/                  # Header files (modular organization)
│   ├── analysis/            # Waveform analysis algorithms
│   ├── config/              # Configuration management
│   └── utils/               # Utility functions
├── src/                     # Implementation files
│   ├── analysis/            # Analysis module implementations
│   ├── utils/               # Utility module implementations
│   └── *.cpp                # Main stage executables
├── output/                  # Output directory (auto-created)
│   ├── root/                # ROOT files
│   ├── hdf5/                # HDF5 files
│   └── waveform_plots/      # Generated plots
├── Makefile                 # Build system
├── converter_config.json    # Unified configuration
└── run_full_pipeline.sh     # Pipeline orchestration script
```

## Pipeline Stages

### Stage 1: waveform_converter (convert_to_root)

Converts binary or ASCII waveform files to ROOT format with pedestal subtraction.

**Input**: Binary files (wave_0.dat, wave_1.dat, ...) or ASCII text files
**Output**: ROOT file in `output/root/` with TTree containing:
- Raw waveforms (ch00_raw, ch01_raw, ...)
- Pedestal-subtracted waveforms (ch00_ped, ch01_ped, ...)
- Time axis, pedestals, board IDs, event counters

**Configuration**: `converter_config.json` → `common` + `waveform_converter` sections
- **Common**: `output_dir`, `n_channels`, `max_cores` (parallel if > 1)
- **waveform_converter**: `input_dir`, `input_pattern`, `input_is_ascii`, `tsample_ns`, `pedestal_window`, `ped_target`

**Key Features**:
- Binary and ASCII format support
- Automatic pedestal subtraction
- Parallel file loading (multi-threaded when `max_cores > 1`)
- Special channel override (e.g., different trigger file)

**Usage**:
```bash
./convert_to_root                    # Uses converter_config.json
./convert_to_root --config my_config.json
./convert_to_root --help
```

### Stage 2: waveform_analyzer (analyze_waveforms)

Extracts timing and amplitude features from waveforms, following the analysis method from ALICE TOF test beam analysis.

**Input**: ROOT file from Stage 1 (in `output/root/`)
**Output**: ROOT file in `output/root/` with analyzed features per channel:

**Extracted Features** (per channel, per event):
- **Baseline & Noise**: baseline, rmsNoise, noise1Point, ampMinBefore, ampMaxBefore
- **Signal**: ampMax, charge, signalOverNoise, peakTime
- **Timing**: riseTime, slewRate
- **CFD Timing**: ch##_timeCFD_##pc, ch##_jitterCFD_##pc (multiple thresholds: 10%, 20%, 30%, 50%)
- **Leading Edge**: ch##_timeLE_##mV, ch##_jitterLE_##mV, ch##_totLE_##mV (multiple voltage thresholds)
- **Charge-based**: ch##_timeCharge_##pc (time to reach % of total charge)

**Configuration**: `converter_config.json` → `common` + `waveform_analyzer` sections
- **Common**: `output_dir`, `n_channels`, `max_cores` (parallel if > 1)
- **waveform_analyzer (per-channel arrays)**: `signal_region_min/max`, `charge_region_min/max`, `baseline_region_min/max`, `signal_polarity`
- **waveform_analyzer (global)**: `cfd_thresholds`, `le_thresholds`, `charge_thresholds`, `impedance`, `sensor_mapping`

**Key Features**:
- Modular analysis functions (in `src/analysis/waveform_math.cpp`)
- CFD and Leading Edge timing with sub-sample precision
- Charge integration and timing
- Per-channel configuration (different regions per channel)
- Parallel processing (multi-process when `max_cores > 1`)
- Optional waveform plotting

**Usage**:
```bash
./analyze_waveforms                  # Uses converter_config.json
./analyze_waveforms --config my_config.json
./analyze_waveforms --help
```

### Stage 3: export_to_hdf5

Exports ROOT data to HDF5 format for analysis in Python/NumPy/pandas.

**Input**: ROOT file from Stage 1 (raw) or Stage 2 (analysis)
**Output**: HDF5 file in `output/hdf5/` for Python/NumPy/pandas analysis

**Two Export Modes**:

1. **Raw Mode** (`--mode raw`):
   - Exports pedestal-subtracted waveforms
   - HDF5 datasets: `Metadata` (event info), `Waveforms` (2D array), `TimeAxis_ns`
   - One row per (event, channel) pair
   - File size: ~33 MB (523 events, 16 channels)

2. **Analysis Mode** (`--mode analysis`):
   - Exports extracted features only (much smaller)
   - HDF5 dataset: `AnalysisFeatures` (structured array)
   - All timing and amplitude features in flat structure
   - File size: ~400 KB (523 events, 16 channels)

**Key Features**:
- Flattened storage for easy Python indexing
- Uses `#pragma pack` for C struct → HDF5 mapping
- Organized output directory structure

**Usage**:
```bash
# Export raw waveforms
./export_to_hdf5 --mode raw

# Export analysis features
./export_to_hdf5 --mode analysis

# Custom files
./export_to_hdf5 --mode analysis --input my_analysis.root --output my_output.hdf5

./export_to_hdf5 --help
```

## Wrapper Script

`run_full_pipeline.sh` runs all stages automatically:

```bash
# Run full pipeline
./run_full_pipeline.sh

# Run specific stages
./run_full_pipeline.sh --stage1-only
./run_full_pipeline.sh --stage2-only
./run_full_pipeline.sh --stage3-only

# Skip stages
./run_full_pipeline.sh --skip-stage2  # Only convert and export, skip analysis

# Custom files
./run_full_pipeline.sh \
  --converter-config my_config.json \
  --analysis-config my_analysis.json \
  --raw-root my_raw.root \
  --analysis-root my_analyzed.root

# Help
./run_full_pipeline.sh --help
```

## Configuration

The pipeline uses a **unified configuration file** `converter_config.json` with hierarchical structure:

### converter_config.json Structure

```json
{
  "common": {
    "output_dir": "output",
    "n_channels": 16,
    "max_cores": 8,              // Auto-enables parallel if > 1
    "chunk_size": 100,
    "temp_dir": "./temp_analysis",
    "waveforms_root": "waveforms.root",
    "waveforms_tree": "Waveforms",
    "analysis_root": "waveforms_analyzed.root",
    "analysis_tree": "Analysis"
  },
  "stage1": {
    "input_pattern": "wave_%d.dat",
    "input_dir": "../AC_LGAD_TEST/",
    "input_is_ascii": false,
    "special_channel_file": "TR_0_0.dat",
    "special_channel_index": 3,
    "enable_special_override": true,
    "tsample_ns": 0.2,
    "pedestal_window": 100,
    "ped_target": 3500.0
  },
  "stage2": {
    "analysis_region_min": [0.0, 0.0, ...],    // Per-channel (16 values)
    "analysis_region_max": [190.0, 190.0, ...],
    "baseline_region_min": [0.0, 0.0, ...],
    "baseline_region_max": [40.0, 40.0, ...],
    "signal_region_min": [0.0, 0.0, ...],
    "signal_region_max": [190.0, 190.0, ...],
    "charge_region_min": [0.0, 0.0, ...],
    "charge_region_max": [190.0, 190.0, ...],
    "signal_polarity": [-1, -1, ...],          // ±1 per channel
    "cfd_thresholds": [10, 20, 30, 50],        // Global (percent)
    "le_thresholds": [10.0, 20.0, 50.0],       // Global (mV)
    "charge_thresholds": [10, 20, 50],         // Global (percent)
    "sensor_mapping": { ... }                  // Maps channels to sensors
  }
}
```

**Key Features**:
- **Hierarchical structure**: `common` section inherited by all stages
- **Per-channel arrays**: Different analysis regions per channel
- **Automatic parallel**: Set `max_cores > 1` to enable parallel processing
- **Fast parsing**: Uses simdjson library for zero-copy JSON parsing
- **Type-safe**: Helper functions in `include/utils/json_utils.h`

## Modular Architecture

The project follows professional C++ organization with clear separation of concerns:

### Source Files (`src/`)
- **Main Executables**:
  - `convert_to_root.cpp` (1146 lines) - Stage 1
  - `analyze_waveforms.cpp` (608 lines) - Stage 2
  - `export_to_hdf5.cpp` (678 lines) - Stage 3
- **Analysis Module** (`src/analysis/`):
  - `waveform_math.cpp` (401 lines) - Core analysis algorithms
  - `waveform_plotting.cpp` (207 lines) - ROOT plotting
- **Utils Module** (`src/utils/`):
  - `file_io.cpp` (185 lines) - Binary/ASCII I/O

### Header Files (`include/`)
- **Config Module** (`include/config/`):
  - `common_config.h` (19 lines) - Base configuration struct
  - `wave_converter_config.h` (119 lines) - Stage 1 config + loader
  - `analysis_config.h` (247 lines) - Stage 2 config + loader
- **Analysis Module** (`include/analysis/`):
  - `waveform_math.h` (108 lines) - Analysis function declarations
  - `waveform_plotting.h` (14 lines) - Plotting interface
- **Utils Module** (`include/utils/`):
  - `json_utils.h` (144 lines) - simdjson wrappers
  - `file_io.h` (42 lines) - File I/O structures
  - `filesystem_utils.h` (68 lines) - Directory/path utilities

### Output Organization (`output/`)
All outputs organized by type:
- `output/root/` - ROOT files (waveforms.root, waveforms_analyzed.root)
- `output/hdf5/` - HDF5 files (waveforms_raw.hdf5, waveforms_analyzed.hdf5)
- `output/waveform_plots/` - Generated plots (if enabled)

**Total**: ~5,600 lines of C++ code across 20 files

## Output Files

| File | Description | Size (typical) |
|------|-------------|----------------|
| `waveforms.root` | Raw waveforms with pedestal subtraction | ~46 MB (523 events) |
| `waveforms_analyzed.root` | Analyzed features | ~250 KB |
| `waveforms_raw.hdf5` | Raw waveforms (HDF5) | ~33 MB |
| `waveforms_analyzed.hdf5` | Analyzed features (HDF5) | ~400 KB |

## Python Analysis Example

```python
import h5py
import numpy as np
import pandas as pd

# Load analysis features
with h5py.File('waveforms_analyzed.hdf5', 'r') as f:
    features = f['AnalysisFeatures'][:]

# Convert to pandas DataFrame
df = pd.DataFrame(features)

# Filter by channel
ch0 = df[df['channel'] == 0]

# Plot amplitude distribution
import matplotlib.pyplot as plt
plt.hist(ch0['ampMax'], bins=50)
plt.xlabel('Amplitude (V)')
plt.ylabel('Counts')
plt.show()

# Analyze timing
print(f"Mean rise time: {ch0['riseTime'].mean()*1e9:.2f} ns")
print(f"Mean S/N: {ch0['signalOverNoise'].mean():.1f}")
```

## Build System

```bash
make               # Build all executables
make clean         # Remove executables
make cleanall      # Remove executables and output files
make test          # Run test pipeline
make help          # Show help
```

## Requirements

- **ROOT 6.x** - CERN data analysis framework
- **HDF5 library** - High-performance data format
  - macOS: `brew install hdf5` → `/opt/homebrew/opt/hdf5`
  - Linux: Usually in `/usr` or `/usr/local`
- **simdjson** - Fast JSON parsing library
  - macOS: `brew install simdjson` → `/opt/homebrew/lib`
  - Linux: Build from source or install via package manager
- **C++17 compiler** - g++ or clang++
- **Platform**: macOS or Linux

## Troubleshooting

**ROOT not found**:
```bash
source /opt/root6/root_install/bin/thisroot.sh
```

**HDF5 not found**:
```bash
# Find HDF5 installation
brew --prefix hdf5  # macOS
find /usr -name "hdf5.h"  # Linux

# Edit Makefile and set HDF5_PREFIX
```

**simdjson not found**:
```bash
# Install simdjson
brew install simdjson  # macOS

# Or edit Makefile INCLUDES and JSON_LIBS paths
```

**Input files not found**:
- Check `input_dir` in `converter_config.json` → `stage1` section
- Ensure path is correct (absolute or relative)
- Verify file naming matches `input_pattern` (e.g., "wave_%d.dat")
