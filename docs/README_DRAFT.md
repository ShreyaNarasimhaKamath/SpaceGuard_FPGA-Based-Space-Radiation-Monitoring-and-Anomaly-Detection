<!-- DRAFT for review. Not committed. Replace README.md only after you have checked every claim and filled the [BRACKETS]. -->
# SpaceGuard: A Self-Protecting Radiation Monitor for FPGA Electronics

SpaceGuard is an ECE engineering project on a radiation-aware detector that runs *inside* an FPGA, watches particle-flux measurements, and decides on its own when the electronics around it should switch to a more protective mode. The detector is designed to protect itself against single-event upsets (SEUs), and its architecture is not tied to one mission: the same structure could be re-tuned for detector electronics at particle accelerators.

The algorithms are developed and evaluated on real GOES-16 particle-flux data, then implemented as board-independent Verilog. The first hardware target is a Xilinx Zynq UltraScale+ (ZCU104).

> **Status.** Everything below is verified in **simulation and Python only**. No design has been run on a physical FPGA, and nothing has been radiation tested. See [Limitations](#limitations).

## What has been built

```text
 counts / flux samples (valid flag)
        |
  per-channel detector  ----  TMR: 3 state copies, whole-state vote every cycle
   (relative-threshold EMA  OR  hazard count accumulator)
        |
  k-of-N channel vote   (are enough independent channels in alarm?)
        |
  mitigation controller  --  NORMAL / ALERT (fast scrub) / PROTECT
        |                     staged de-escalation, watchdog fail-safe on missing data
  scrub_fast, protect   (outputs that drive the FPGA's own protection)
```

- **Two detector types** (`rtl/anomaly_detector/`): a shift-only asymmetric exponential baseline with a relative threshold (for fast-varying flux such as the 50-120 keV electron bands), and a leaky count accumulator with an absolute alarm level (for hazard channels such as the >2 MeV electron integral flux).
- **Channel vote and controller** (`rtl/top/`): `sg_top.v` combines N channels with a k-of-N vote; `sg_ctrl.v` turns the vote into a protection level.
- **Self-protection:** detector and controller state are triplicated with a bitwise majority vote whose result is fed back to all copies, so a single upset is overwritten on the next clock edge.
- **Golden models** (`python/`): bit-exact integer models of the detectors and controller generate the RTL test vectors.

## Results so far

All detection results use synthetic injected events on real data unless stated otherwise; the three days of data contain no confirmed hazard event.

| Item | Result |
|---|---|
| Relative-threshold EMA vs the earlier linear-threshold detector (held-out day, injected events) | detection of core events 1.00 vs 0.96, natural alarms/day 3.5 vs 4.0 (Python study, `docs/v2_findings.md`) |
| Fixed-point implementation | 8 fractional bits reproduce the floating-point alarms exactly on this data |
| RTL vs golden model | 0 alarm, baseline or mode mismatches over 86,400-sample vectors (clean, injected, dropout/blackout) |
| Single-bit upsets in detector/controller state (169 per run, 3 seeds) | unprotected: 20-54 events caused a wrong output; TMR version: 0. Same bit hit in two copies defeats TMR (flagged but not corrected) |
| Hazard channel, real data | >2 MeV integral flux never reaches 1000 pfu on 5-minute averages; 0 alarms at a 1000 pfu alert level on 6 channel-days |
| Hazard channel, injected Poisson bursts | all bursts at 2x the 1000 pfu level or above detected, median delay about 150 s at 2x and 48 s at 10x; marginal events (about 1000 pfu for 5 minutes) often missed |
| Resource use, Zynq UltraScale+ (xczu7ev, Vivado 2025.2, out-of-context, 300 MHz target) | relative-threshold EMA system (2 channels + vote + controller): bare 632 LUT / 162 FF, TMR 2211 LUT / 486 FF (0.96% of LUTs). Count (hazard) system: bare 257 LUT / 96 FF, TMR 959 LUT / 288 FF (0.42% of LUTs). No BRAM or DSP. TMR costs about 3.5-3.7x LUTs and exactly 3x flip-flops. Timing at 300 MHz: bare variants meet it; TMR EMA misses by 0.05 ns (about 295 MHz), TMR count meets it (about 400 MHz). Not a board-level result: no I/O, no place-and-route with real constraints |

## Data

**GOES-16 SEISS MPS-HI Level 1b, daily aggregations** (NOAA). Telescopes 1 (and 4 for the hazard channel), 1 Hz, **1-3 January 2024**.

> NOAA National Centers for Environmental Information. *GOES-R Series Level 1b SEISS Magnetospheric Particle Sensor High Energy (MPS-HI) Product in Daily Aggregations* (GOES-16). [Add: URL of the exact download page you used, and access date: DD Month YYYY].
> Product page: https://ncei.noaa.gov/products/goes-r-space-environment-in-situ

Raw `.nc` files are not stored in this repository. Samples are used according to the file's data-quality flags (details in `docs/dataset_audit.md`). Results here are the author's own analysis and are not endorsed by NOAA or NASA.

## Reproducing the results

```text
pip install numpy pandas scipy netCDF4 numba matplotlib
cd python
python dataset_audit.py                      # data checks
python v2_injection_benchmark.py             # EMA detector benchmark
python v2_fixedpoint_upset.py                # fixed-point + upset study
python v2_hazard_counts.py                   # hazard channel study
python v2_rtl_refmodel.py ../results/rtl_vectors      # RTL test vectors
python v2_system_model.py ../results/rtl_vectors
python v2_hazard_vectors.py ../results/rtl_vectors
```

RTL simulation uses Icarus Verilog, for example:
`iverilog -Ptb_top.K=2 -o sim rtl/anomaly_detector/*.v rtl/top/*.v testbench/tb_top.v && vvp sim +PFX=results/rtl_vectors/sys_inj`
Vivado synthesis (out of context, ZCU104 part): `vivado\run_ooc.bat` or `source vivado/run_ooc.tcl` in the Vivado Tcl console.

## Limitations

- Simulation only: no hardware run, no radiation or fault-injection testing on a device.
- Upsets modelled: single bit flips in the detector and controller state registers. **Not modelled:** multi-bit upsets, upsets in the voter, transients in combinational logic, configuration-memory upsets (the dominant effect in SRAM FPGAs), latch-up.
- Three days of quiet data, one hazard channel type; detection performance uses synthetic events. The MPS-HI data has no >=10 MeV proton channel, so proton-storm hazards are not covered.
- Counts for the hazard channel are reconstructed from flux (about 1326 pfu per count/s), so they are a model of the raw counter data.
- Detection thresholds and controller timings (ESC, HOLD, WD) are design parameters, not mission-validated values.

## Repository layout

```text
docs/       dataset audit, v2 findings, RTL notes, this README draft
python/     analysis, golden models, vector generators, studies
rtl/        anomaly_detector/ (detectors), top/ (channel vote + controller)
testbench/  exactness and fault-injection testbenches
vivado/     out-of-context synthesis script for the ZCU104 part
results/    generated outputs (not all tracked)
```

## Roadmap

- [x] Detector design, fixed-point model, RTL, TMR, fault-injection simulation
- [ ] Vivado resource/timing/power numbers (bare vs TMR)
- [ ] Decide the protective action driven by `protect` on the ZCU104 (for example the Soft Error Mitigation controller)
- [ ] Further upset models (multi-bit, voter, state-machine corruption)
- [ ] Board-level test on the ZCU104

## Acknowledgements

- **Data:** NOAA National Centers for Environmental Information and the NOAA GOES-R program for the freely available GOES-16 SEISS MPS-HI data.
- **Tools:** Python (NumPy, pandas, SciPy, netCDF4, Numba, Matplotlib), Icarus Verilog, and AMD Vivado.
- **Guidance:** [Faculty guide / mentor name, department], VIT Chennai. [Add or delete]
- **AI assistance:** [Optional, check your department's policy: "Parts of the analysis scripts, RTL and documentation were developed with the help of an AI assistant (Anthropic's Claude); all results were reviewed by the author."]

## Author

[Your name], Electronics and Communication Engineering, VIT Chennai.

## License

[Not chosen yet. Pick one (for example MIT or Apache-2.0) before publishing, or state "all rights reserved".]

## Disclaimer

An academic engineering project. Results come from historical satellite data and simulation and must not be used as operational radiation alerts or as evidence of flight-qualified fault tolerance.
