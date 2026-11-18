# SpaceGuard dataset audit — GOES-16 SEISS MPS-HI L1b, 2024-01-01 to 2024-01-03

Reproduce with `python python/dataset_audit.py` (outputs in `results/dataset_audit/`). Not committed to git.

## 1. What the files are
- Daily NCEI aggregation of 30 s L1b files (2,880 per day); 86,400 reports/day, timestamps are the *end* of each 1 s observation.
- 5 telescopes x 10 differential electron bands (+ >2 MeV integral channel E11), 5 x 11 proton bands, 2 dosimeters (LoLET, HiLET), plus eclipse / mode / processing flags and spacecraft position.
- Band labels in the file are generic ("ElectronEnergyBand-N"). NOAA describes the 10 bands as log-spaced over 50 keV-4 MeV. If exactly log-spaced, Band 1 ~ 50-78 keV and Band 2 ~ 78-120 keV (derived here, not read from the file; confirm against the product user guide).

## 2. Timing and instrument state
- Perfect 1 Hz on all three days: no gaps, duplicates or non-1 s steps.
- Jan 2 has 178 consecutive missing rows (all variables fill / mode flag missing / processing flag 5).
- Instrument always in operational mode; eclipse flag always 0 (January is outside eclipse season, so eclipse effects are untested).

## 3. Quality flags (DQF is a bitmask, not a boolean)
Values seen: 0 good, 8 out-of-band contamination, 16 dynamic-error threshold exceeded, 24 both, 255 fill.
- DQF 8 <=> negative flux (859,902 samples; no other cause). DQF 16 <=> small *positive* flux.
- DQF 0 contains 2.9 M exact-zero fluxes. So `dqf == 0` keeps zeros and discards small positive values: a low-flux-biased gap in the record.
- Bands 1-2: ~99-100 % good on all telescopes except Band 2 on T2/T3/T5 (several % flagged 16). On Jan 2 the flagged Band 2 samples sit at lower flux than the surrounding good data (local median ~780 vs day median 6,869).
- Bands 3-10: validity is 1-100 % and varies strongly by day/band; whether a sample is "valid" is tied to the flux level. Do not use `dqf == 0` masks on these bands without handling this.
- >2 MeV integral channel: valid samples are all zero on all three days.

## 4. Noise and variability
- Reported uncertainty is ~19-21 % of flux, almost constant, while the real 1 s scatter is only 0.11-0.27 of it. It is not a counting-noise estimate and cannot be used as a per-sample noise model.
- Short-term noise grows ~ sqrt(level): ~5 % of level at low flux, ~1.4 % at high flux (Band 1, Jan 1). A fixed absolute threshold is mismatched to the signal.
- std(log10 flux) of block means is almost unchanged from 1 s to 30 min (e.g. 0.427 -> 0.414). Variability is slow, large-scale and non-stationary; 1 s noise is negligible for Bands 1-2.
- The five telescopes track each other (Band 1, 1-min means: median pairwise correlation 0.84-0.96), so the large changes are coherent, not single-channel glitches. NOAA documents telescope-to-telescope scale factors of ~0.68-2.3 in its read-me, so absolute levels across telescopes are not interchangeable.

## 5. Residual stationarity (post-10 h MAD / first-10 h MAD, rolling-median baseline)
| day / band | linear | relative (x/B-1) | log10 |
|---|---|---|---|
| Jan 1 B1 | 12.58 | 3.02 | 3.06 |
| Jan 1 B2 | 9.77 | 2.38 | 2.47 |
| Jan 2 B1 | 3.36 | 0.80 | 0.82 |
| Jan 2 B2 | 2.79 | 0.73 | 0.74 |
| Jan 3 B1 | 2.92 | 2.04 | 2.03 |
| Jan 3 B2 | 1.77 | 1.86 | 1.82 |
Relative/log residuals are markedly more stable in 4 of 6 cases and about equal in the other two.

## 6. Radiation-hazard relevance of the development days
- Bands 1-2 (~50-120 keV electrons) are a magnetospheric-dynamics signal; they are not the quantities used for SEU / deep-dielectric-charging alerts.
- Operational definitions use >2 MeV electron flux > 1000 pfu and >=10 MeV proton flux >= 10 pfu (NOAA SWPC). Valid >2 MeV electron flux is zero on all 3 days; the highest proton bands have zero median. Integral >=10 MeV proton flux was not computed here. The three days look quiet for hazard channels — to be confirmed against the SWPC alert archive.
