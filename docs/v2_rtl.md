# SpaceGuard v2 - single-channel RTL and fault injection

Files: `rtl/anomaly_detector/sg_next.v` (next-state logic), `sg_detector.v` (simplex and TMR tops),
`testbench/tb_exact.v`, `testbench/tb_fault.v`, `python/v2_rtl_refmodel.py` (integer golden model + vector generator),
vectors in `results/rtl_vectors/`, log in `results/v2/rtl_fault_injection.log`.
Simulator: Icarus Verilog 11.0 (no vendor tool). Plain Verilog-2001, board independent.

## Design
66-bit state: baseline 32, freeze counter 11, off 6, run 4, warm counter 11, on, init. Shift-only EMA (up 2^-10, down 2^-7),
threshold base + base>>1, M=10, R=30, warm-up 1800, freeze <=1800. TMR: three copies, bitwise vote of the whole state every cycle,
voted state fed back to all copies (an upset is overwritten the next edge). `tmr_err` flags any copy that disagrees.

## Bit-exact check (RTL vs integer model, which also equals the numba study model)
3 vector sets x 86,400 samples (injected Jan 3 B1, clean Jan 1 B2, injected Jan 2 B1 with 60 dropout bursts):
0 alarm mismatches and 0 baseline mismatches, simplex and TMR.

## Fault injection (one state-bit upset every 500 samples, 169 events per run, 3 seeds, 2 vector sets)
- Single upset in one copy: simplex had wrong alarm in 29-48 of 169 events (about 10-16k wrong cycles); TMR 0 events, 0 wrong cycles; tmr_err flagged all 169.
- Same bit hit in two copies at once: TMR identical to simplex (voter defeated), but tmr_err still flags it.

## Not covered / caveats
- Upsets are injected on clock-edge boundaries into state flops only; no combinational transients (SET), clock/reset upsets, config memory, or multi-bit upsets across different bits.
- The voter and tmr_err logic are not themselves faulted. Simplex/TMR LUT/timing numbers are not measured (Yosys LUT mapping not completed); TMR is 198 state flops by construction.
- Upset rate (every 500 samples) is stress testing, not a physical rate.
- Nothing committed to git.

## System level: channel vote + mitigation controller (added)
Files: `rtl/top/sg_top.v`, `rtl/top/sg_ctrl.v`, `testbench/tb_top.v`, `tb_top_fault.v`, `python/v2_system_model.py`; log `results/v2/rtl_system_fault_injection.log`.
- N=2 channels (electron B1, B2) -> k-of-N vote (K=1 or 2) -> controller. Modes NORMAL / ALERT (scrub_fast) / PROTECT (protect).
  ESC=300 ticks of event -> PROTECT; staged de-escalation with HOLD=600 ticks per stage; watchdog WD=60 ticks with no valid sample on any channel -> PROTECT (fail-safe). Illegal mode 3 behaves as PROTECT. HOLD/ESC/WD are placeholders, not tuned.
- Controller state (30 bits) is TMR with whole-state vote and feedback, same scheme as the detectors.
- Bit-exact vs `v2_system_model.py` for K=1 and K=2 on an injected day and a dropout/blackout day: 0 mismatches.
- Fault injection (one upset per 500 samples into a random detector channel or the controller, 169 events, 3 seeds, 2 sets): bare top had wrong mode in 20-43 events (6.9k-14.6k wrong cycles); protected top 0. Same bit in two copies: protected top = bare top.
- Controller sees detector alarms with 1-sample latency.

### Caveat that matters: occupancy on CLEAN data
Mode occupancy of the un-injected days (samples of 86,400 in NORMAL/ALERT/PROTECT): Jan1 K=2 72677/4500/9223, Jan2 K=2 73715/2836/9849, Jan3 K=2 80252/1800/4348 (K=1 is higher: PROTECT 11k-19k). The detector alarms on natural variation of 50-120 keV electrons, so with these placeholder timings the controller would sit in PROTECT 5-11% of the day. The mitigation policy therefore needs hazard-relevant channels (higher-energy electrons, protons) and tuned ESC/HOLD, not B1-B2 alone.

## Hazard channels: >2 MeV electron integral as a count detector (added)
Files: `python/v2_hazard_counts.py`, `v2_hazard_vectors.py`, `rtl/anomaly_detector/sg_cdet.v`, `sg_top.v` (KIND=1), `testbench/tb_top_hz.v`, `tb_top_fault_hz.v`;
results `results/v2/hazard_*.csv`, `rtl_hazard_fault_injection.log`.
- Channel: IntgElectronFluxes, telescopes 1 and 4 (both documented as agreeing for >2 MeV). Input to the RTL is integer counts/s, c = round(flux/1326); negative (background-subtracted) samples count as 0 and stay valid; only NaN/fill is invalid. Quantisation is approximate (residual up to 0.5 count), so this is a model of the counts, not raw counter data.
- The 3 days never reach the SWPC >2 MeV alert level: max 5-min mean 323-433 pfu (max 433), 0 of 288 bins >= 1000 pfu on every channel-day. Mean background about 0.17 counts/s (about 230 pfu).
- Detector: L <- L - (L>>8) + (c<<8) with 8 fractional bits, absolute Schmitt alarm (TH 1000 pfu = 49424, TL 70% = 34597). Without fractional bits the floor shift stops decay below 256 (found and fixed: it caused false alarms).
- Real data alarms by alert level: 0 at 500-2000 pfu; 2 at 400 pfu.
- Poisson injection on real counts (24 events per cell): added rate A x 1000 pfu on top of background. Detected 100% of A>=2 events at all durations, median delay 150 s at 2x, 72 s at 5x, 48 s at 10x; marginal (A=1, 300 s) only 38%/8% (1of1/2of2); A<=0.5 (total below 1000 pfu) no alarms. A=0.8 for 300 s is also missed. 2-of-2 vote adds roughly 20-100 s delay vs 1-of-1.
- System (2 count channels -> vote -> controller): clean real day stays NORMAL for all 86,400 samples; RTL bit-exact vs model (K=1,2; clean and injected+dropout sets); fault injection same as before (bare 41-54 of 169 events wrong, protected 0; same-bit double hit defeats TMR).
- Limits: only one hazard channel type; no proton >=10 MeV channel exists in MPS-HI (protons span 80 keV-12 MeV; the top 3 channels register only in solar particle events, per NOAA readme), so SEU-relevant proton hazards need SGPS/EHIS data that we do not have. Injection ground truth is synthetic. Alert level is a parameter, not tuned for a specific mission.

## Controller timing study (added)
`python/v2_ctrl_tuning.py`, `results/v2/controller_timing_study.csv`. Hazard system, 2-of-2 vote, 15 hazard-level injected events and 3 sub-hazard events (3 days).
- All 15 hazard events reach ALERT at every setting (median 168 s after event start, dominated by the detector's ~150 s response, not by ESC). ESC sets only the PROTECT latency: about +60/+120/+300/+600 s after ALERT for ESC=60/120/300/600; at ESC=600 three of 15 events never reach PROTECT (shorter than ESC).
- HOLD sets the hangover: median samples non-NORMAL after an event ends = about 2 x HOLD + 320 s (560 s at HOLD=120, 1520 s at 600, 3920 s at 1800), because de-escalation goes PROTECT -> ALERT -> NORMAL with HOLD each. No flapping at any setting (4 mode changes per event).
- Non-NORMAL samples on the clean days (358 at HOLD=120 up to 3718 at 1800) come from one real data gap, Jan 2 samples 61928-62105 where both telescopes are invalid, which triggers the watchdog fail-safe; they are not radiation alarms. The one 'sub-hazard trigger' at HOLD=1800 is the hangover from that same gap reaching the sub-hazard event window (62105 + 2 x 1800 > 64000).
- Policy choice, not an optimum: ESC about 300 s matches the 5-minute averaging SWPC uses; HOLD 300-600 s keeps the hangover to 15-25 minutes. Current defaults (ESC 300, HOLD 600) are inside that range. Final values depend on what PROTECT costs the system.

## Vivado 2025.2 out-of-context results, xczu7ev-ffvc1156-2-e (ZCU104/106), 300 MHz target (added)
Script `vivado/run_ooc.tcl`; reports and `summary_xczu7ev-ffvc1156-2-e.csv` in `results/vivado/`. Full system (2 channels + vote + controller).
| variant | LUT | FF | CARRY8 | WNS at 300 MHz (ns) | est. fmax (MHz) | total / dynamic power (W) |
|---|---|---|---|---|---|---|
| bare_ema | 632 | 162 | 46 | +0.220 | 321 | 0.601 / 0.010 |
| tmr_ema | 2211 | 486 | 162 | -0.051 | 296 | 0.636 / 0.044 |
| bare_count | 257 | 96 | 10 | +1.459 | 534 | 0.598 / 0.006 |
| tmr_count | 959 | 288 | 48 | +0.844 | 402 | 0.616 / 0.024 |
- TMR = 3x flip-flops exactly (486 = 2 detectors x 3 x 66 + 3 x 30; 288 = 2 x 3 x 33 + 90), so the DONT_TOUCH attributes kept all copies. LUTs 3.5x (EMA) and 3.7x (count). Under 1% of the part; no BRAM/DSP.
- Power: static 0.592 W dominates; dynamic figures are vectorless estimates with default activity, not meaningful for comparison beyond 'small'.
- Caveats: out-of-context (no I/O, no real clock tree/board constraints); timing misses of 0.05 ns at 300 MHz are irrelevant at the 1 Hz sample rate but matter if reused at accelerator rates. A first run (before the controller fix) counted two unused controller copies in the bare variants; those numbers were discarded.
- Not yet checked: that the TMR voters survive in the netlist as separate majority gates, and any effect of the shared (not triplicated) next-state input nets.
