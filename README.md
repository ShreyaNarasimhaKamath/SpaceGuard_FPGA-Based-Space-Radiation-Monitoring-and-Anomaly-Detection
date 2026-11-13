# SpaceGuard: FPGA-Based Space Radiation Monitoring and Anomaly Detection

## Overview

SpaceGuard is an electronics engineering project focused on developing a digital signal processing and FPGA-oriented architecture for spacecraft radiation monitoring and anomaly detection using publicly available NASA datasets.

The project investigates the processing of spacecraft radiation measurements, identification of abnormal radiation patterns, and hardware-efficient implementation of digital monitoring algorithms.

## Objectives

- Analyse real spacecraft radiation measurements from NASA datasets.
- Develop digital signal processing techniques for radiation measurement analysis.
- Investigate statistical and adaptive anomaly detection algorithms.
- Design fixed-point digital processing architectures using Verilog HDL.
- Explore fault-tolerant hardware techniques, including Triple Modular Redundancy (TMR).
- Evaluate detection performance and hardware implementation metrics.

## Proposed System Architecture

```text
NASA Spacecraft Radiation Dataset
              |
              v
     Data Acquisition
              |
              v
    Signal Preprocessing
              |
              v
       Digital Filtering
       (FIR / IIR)
              |
              v
     Anomaly Detection
   (Z-score / CUSUM)
              |
              v
    FPGA Hardware Design
              |
              v
     Radiation Event Alert
```

## Dataset

Initial dataset: NASA HEO-1 EXP radiation measurements.

The dataset contains particle-count and radiation-dose-related measurements collected by a spacecraft instrument.

Dataset source:

https://data.nasa.gov/dataset/heo-1-data-from-the-exp-instrument-in-15-second-averages-54cf8

## Technologies

- Python
- NumPy
- Pandas
- SciPy
- Matplotlib
- Verilog HDL
- AMD Vivado (planned)

## Development Roadmap

- [ ] NASA dataset acquisition and exploration
- [ ] Data cleaning and preprocessing
- [ ] Statistical characterisation of radiation measurements
- [ ] Digital filtering and signal analysis
- [ ] Radiation anomaly detection
- [ ] Algorithm performance evaluation
- [ ] Fixed-point hardware architecture
- [ ] Verilog RTL implementation
- [ ] RTL simulation and verification
- [ ] FPGA synthesis and resource analysis
- [ ] TMR architecture and fault-injection evaluation

## Project Status

Initial development phase: repository setup and dataset exploration.

## Author

Electronics and Communication Engineering

## Disclaimer

This is an academic engineering project. Radiation detection and fault-tolerance capabilities will be evaluated through documented software experiments, RTL simulation and physical hardware testing where available.