# SpaceGuard v2 - injection benchmark, fixed-point and upset study (working notes, not committed)

Scripts: `python/v2_injection_benchmark.py`, `python/v2_fixedpoint_upset.py`. Outputs: `results/v2/`.
Data: GOES-16 MPS-HI, Telescope 1, electron Bands 1-2, 2024-01-01..03 only. Dev = Jan 1-2, held-out test = Jan 3.

## Method
- Ground truth: trapezoidal enhancements (amplitude 0.5/1/2/4 x added, plateau 60/300/1800 s, rise 10/120 s) multiplied onto the real flux; 8 events/day/rep, 25 reps per series, random positions (different seeds for dev and test).
- "Natural alarms/day" = alarms on the un-injected record. These are NOT counted as false alarms, because the real record contains genuine large excursions. The trade-off curve (detection probability vs natural alarms/day) is the figure of merit.
- Masking: keep DQF 0 and 16 (low-count degraded), drop 8, 24, 255.
- Detectors are causal: median(1800) / asymmetric EMA (alpha = 2^-n) with relative threshold x > (1+k)*baseline; persistence 10 s, release 30 s; optional baseline freeze during alarm (<=1800 s). Legacy = old linear 6-sigma median detector, same persistence.

## Benchmark (core = amplitude >=2x, duration >=5 min)
| config | dev Pd core | dev nat. alarms/day | test Pd core | test nat. alarms/day |
|---|---|---|---|---|
| LEGACY linear 6 sigma | 0.89 | 8.3 | 0.96 | 4.0 |
| MED k=0.5 | 0.94 | 6.5 | 0.97 | 6.5 |
| EMA 2^-10/2^-7, k=0.5, freeze | 1.00 | 6.3 | 1.00 | 3.5 |
| EMA 2^-10/2^-7, k=0.75, no freeze | 0.96 | 5.0 | 0.99 | 2.0 |
Median alarm delay 17-24 s (dominated by the 10 s persistence + ramp). Freeze raises alarm coverage of an event from ~0.7 to ~0.9 and raises natural alarm duty cycle.

## Fixed point (32-bit baseline register, alpha shifts 10/7, k=0.5 as s+(s>>1))
Fractional bits 0/2/4/8/12 -> alarm agreement with float 0.940 / 0.985 / 0.9994 / 1.0 / 1.0. 8 fractional bits are bit-exact on this data. Flux needs ~19 integer bits.

## Single-event upsets of the detector's own state (one bit flip, errors counted over 3 h)
- EMA baseline register: flips with weight below ~2^2 flux units are harmless (self-healing); flips with weight >= 2^10 cause an alarm-state error in ~20-37 % of cases, p90 ~ 550-730 s (no freeze) and ~2,200 s (freeze: the corrupted baseline is held frozen).
- Median buffer: flipping one stored sample has essentially no effect (max 1 s), but the exposure is 1800 words x 32 b ~ 57,600 bits vs ~32 + ~60 control bits for the EMA path.
- Not modelled: multi-bit upsets, upsets in counters/FSM/comparators, configuration-memory upsets (dominant in SRAM FPGAs), SEFI/latch-up.

## Caveats
Idealised injected shapes; one telescope; Bands 1-2 only; 3 days; k, M, R, tmax not jointly optimised; test day is a single day; no accelerator-rate data.
