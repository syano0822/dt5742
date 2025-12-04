# CAEN DT5742 Waveform Processing Pipeline

A three-stage pipeline for converting and analyzing CAEN DT5742 digitizer data.

## Quick Start

```bash
# 1. Setup ROOT environment
source /opt/root6/root_install/bin/thisroot.sh # example!

# 2. Build
make

# 3. Edit configuration
# Edit converter_config.json with your input paths and settings

# 4. Run
./run_full_pipeline.sh
```

## Pipeline Overview

```
Binary/ASCII → [Stage 1] → ROOT → [Stage 2] → ROOT + Features → [Stage 3] → HDF5
                convert_to_root    analyze_waveforms              export_to_hdf5
```

**Stage 1**: Converts binary/ASCII waveform files to ROOT with pedestal subtraction
**Stage 2**: Extracts timing and amplitude features (CFD, leading edge, charge, etc.)
**Stage 3**: Exports to HDF5 format for Python/NumPy/pandas analysis

## Configuration

Edit `converter_config.json`:

```json
{
  "common": {
    "output_dir": "output",
    "n_channels": 16,
    "max_cores": 8              // Set > 1 to enable parallel processing
  },
  "stage1": {
    "input_dir": "../AC_LGAD_TEST/",
    "input_pattern": "wave_%d.dat",
    "input_is_ascii": false
  },
  "stage2": {
    "signal_region_min": [0.0, 0.0, ...],
    "signal_region_max": [190.0, 190.0, ...]
    // ... per-channel analysis settings
  }
}
```

## Usage

### Run Full Pipeline
```bash
./run_full_pipeline.sh                # All stages
./run_full_pipeline.sh --stage1-only  # Convert only
./run_full_pipeline.sh --stage2-only  # Analyze only
```

### Run Individual Stages
```bash
./convert_to_root                     # Stage 1
./analyze_waveforms                   # Stage 2
./export_to_hdf5 --mode analysis      # Stage 3
```

## Output Files

All outputs in `output/` directory:

- `output/root/waveforms.root` - Raw waveforms (~46 MB for 523 events)
- `output/root/waveforms_analyzed.root` - Extracted features (~250 KB)
- `output/hdf5/waveforms_analyzed.hdf5` - Features in HDF5 format (~400 KB)

## Python Analysis

```python
import h5py
import pandas as pd

# Load data
with h5py.File('output/hdf5/waveforms_analyzed.hdf5', 'r') as f:
    data = f['AnalysisFeatures'][:]

# Analyze
df = pd.DataFrame(data)
ch0 = df[df['channel'] == 0]
print(f"Mean amplitude: {ch0['ampMax'].mean():.3f} V")
```

## Requirements

- ROOT 6.x
- HDF5 library (`brew install hdf5`)
- simdjson library (`brew install simdjson`)
- C++17 compiler

## Troubleshooting

**ROOT not found**: `source /opt/root6/root_install/bin/thisroot.sh`
**HDF5 not found**: Edit `HDF5_PREFIX` in `Makefile`
**Input files not found**: Check `input_dir` in `converter_config.json`

## Documentation

- **CLAUDE.md** - Detailed documentation (architecture, algorithms, development guide)
- **AGENTS.md** - Agent-based workflow automation
- Each executable has `--help` option

## Project Structure

```
data_converter/
├── include/          # Headers (analysis, config, utils modules)
├── src/              # Implementation files
├── output/           # Generated outputs (root/, hdf5/, plots/)
├── Makefile          # Build system
└── converter_config.json  # Unified configuration
```

Total: ~5,600 lines of modular C++ code across 20 files.
