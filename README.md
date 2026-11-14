# SpaceGuard: FPGA-Based Space Radiation Monitoring and Anomaly Detection

## Overview

SpaceGuard is an ECE engineering project exploring signal processing and FPGA-oriented anomaly detection using publicly available GOES-16 satellite radiation measurements.

The project analyses particle flux measurements from the GOES-16 Space Environment In-Situ Suite (SEIS), develops statistical methods to identify candidate flux anomalies, and investigates how the processing pipeline could be mapped to digital hardware.

**Current implementation status:** Python-based data analysis and detector evaluation. Fixed-point modelling and Verilog RTL implementation are planned; physical FPGA implementation and radiation fault-tolerance validation have not yet been completed.

## Objectives

- Analyse real satellite radiation measurements.
- Validate and preprocess time-series data using instrument quality flags.
- Develop baseline estimation and residual-based anomaly detection.
- Evaluate detection sensitivity across multiple channels and days.
- Investigate hardware-efficient filtering and detection architectures.
- Develop a fixed-point model and Verilog RTL implementation.
- Explore fault-tolerant techniques such as Triple Modular Redundancy (TMR) as a future extension.

## Dataset

**Primary dataset:** GOES-16 SEIS L1b MPS-HI Science Data.

- Source: NOAA National Centers for Environmental Information (NCEI).
- Sampling frequency: 1 Hz.
- Measurement: Differential electron flux.
- Data format: NetCDF (`.nc`).
- Analysis period currently used: 1–3 January 2024.
- Selected channels: Telescope 1, Electron Bands 1 and 2.

Flux measurements are filtered using the associated data-quality flags. Samples are treated as valid only when the quality flag is zero and the measurement is finite and non-negative.

Raw satellite data are not stored in this repository. Refer to `docs/dataset_readme.txt` for acquisition details.

## Processing Architecture

```text
GOES-16 SEIS L1b MPS-HI Data
              |
              v
      Data Validation
              |
              v
     Signal Preprocessing
              |
              v
       Baseline Estimation
              |
              v
      Residual Calculation
              |
              v
   Robust Threshold Detection
              |
              v
       Event Grouping
              |
              v
  Cross-Channel Event Association
              |
              v
   FPGA-Oriented Architecture
              |
              v
       Verilog RTL (Planned)
```

## Current Detection Method

The current reference detector uses:

- A 30-minute rolling-median baseline.
- Residual calculation between measured flux and estimated baseline.
- Robust noise estimation using the median absolute deviation (MAD).
- A provisional threshold of six times the robust sigma estimate.
- A 120-second event-merging gap.
- A minimum of 60 active threshold-crossing samples per event.
- Independent detection on Electron Bands 1 and 2, followed by temporal event association.

These parameters are experimental and are being evaluated for sensitivity and stability. Detected intervals are treated as candidate flux anomalies, not independently confirmed physical radiation events.

## Repository Structure

```text
SpaceGuard/
├── docs/
│   └── dataset_readme.txt
├── data/
│   ├── raw/
│   └── processed/
├── python/
│   └── Analysis and evaluation scripts
├── results/
│   ├── plots/
│   └── evaluation/
├── rtl/
│   ├── filters/
│   ├── anomaly_detector/
│   └── top/
├── testbench/
├── reports/
├── README.md
└── requirements.txt
```

## Technologies

- Python 3.10
- NumPy
- Pandas
- SciPy
- Matplotlib
- netCDF4
- Verilog HDL (planned)
- AMD Vivado (planned)

## Current Progress

- [x] GOES-16 dataset acquisition for three days.
- [x] NetCDF data inspection and channel selection.
- [x] Data-quality validation.
- [x] Signal characterization and visualization.
- [x] Rolling-median residual-based detector.
- [x] Independent channel detection.
- [x] Multiday event comparison.
- [x] Threshold sensitivity analysis.
- [x] Baseline estimator comparison.
- [ ] Finalize a hardware-oriented baseline and detection algorithm.
- [ ] Develop fixed-point reference model.
- [ ] Implement Verilog RTL.
- [ ] Create RTL testbenches and verify functional equivalence.
- [ ] Evaluate synthesis resources and latency.
- [ ] Investigate TMR and fault-injection validation.

## Scientific Limitations

The dataset contains satellite particle-flux measurements; it does not provide ground-truth labels for the anomalies identified by this project. Candidate detections therefore require independent scientific validation before being interpreted as confirmed radiation events.

The current work is software-based. It does not demonstrate actual FPGA radiation tolerance or immunity to single-event upsets.

## Author

Electronics and Communication Engineering  
VIT Chennai

## Disclaimer

This is an academic engineering project. Results represent experimental analysis of historical satellite measurements and should not be interpreted as validated operational radiation alerts.
