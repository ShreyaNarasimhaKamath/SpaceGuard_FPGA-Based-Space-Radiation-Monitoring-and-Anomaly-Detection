"""
SpaceGuard dataset audit (GOES-16 SEISS MPS-HI L1b, 2024-01-01..03).
Reproduces the numbers quoted in docs/dataset_audit.md. Writes CSVs to results/dataset_audit/.
Read-only on the data; does not touch git.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "raw" / "goes16_mpshi"
OUT = ROOT / "results" / "dataset_audit"
OUT.mkdir(parents=True, exist_ok=True)
DATES = ["20240101", "20240102", "20240103"]
CAL = 10 * 3600


def load(d):
    with Dataset(DATA / f"goes16_mpshi_{d}.nc") as nc:
        g = lambda v: np.ma.filled(nc[v][:].astype(float), np.nan)
        gi = lambda v: np.ma.filled(nc[v][:].astype(float), 255).astype(int)
        return dict(t=g("L1a_SciData_TimeStamp"), f=g("DiffElectronFluxes"), u=g("DiffElectronUncertainties"),
                    q=gi("DiffElectronFluxDQFs"), p=g("DiffProtonFluxes"), pq=gi("DiffProtonFluxDQFs"),
                    ie=g("IntgElectronFluxes"), ieq=gi("IntgElectronFluxDQFs"),
                    lo=g("Dos1_LoLetDose"), hi=g("Dos1_HiLetDose"),
                    mode=gi("Instrument_Mode"), ecl=gi("eclipse_flag"), proc=gi("L1b_Processing_Flag"))


def mad(a):
    return np.median(np.abs(a - np.median(a)))


data = {d: load(d) for d in DATES}

# 1. timing and instrument state
rows = []
for d, x in data.items():
    dt = np.diff(x["t"])
    rows.append(dict(date=d, samples=len(x["t"]), non_1s_steps=int((np.abs(dt - 1) > 0.01).sum()),
                     missing_rows=int((x["mode"] == 255).sum()), operational_mode_rows=int((x["mode"] == 2).sum()),
                     eclipse_rows=int(((x["ecl"] != 0) & (x["ecl"] != 255)).sum())))
pd.DataFrame(rows).to_csv(OUT / "timing_and_state.csv", index=False)

# 2. DQF semantics: flag value vs sign of flux, all 50 channels
F = np.concatenate([x["f"] for x in data.values()]); Q = np.concatenate([x["q"] for x in data.values()])
sign = np.where(np.isnan(F), "nan", np.where(F < 0, "neg", np.where(F == 0, "zero", "pos")))
pd.crosstab(Q.ravel(), sign.ravel()).to_csv(OUT / "dqf_vs_flux_sign.csv")

# 3. percent good (dqf==0) per day x telescope x band
rows = []
for d, x in data.items():
    good = (x["q"] == 0) & np.isfinite(x["f"])
    for t in range(5):
        for b in range(10):
            rows.append(dict(date=d, telescope=t + 1, band=b + 1, pct_good=100 * good[:, t, b].mean()))
pd.DataFrame(rows).to_csv(OUT / "pct_good_by_channel.csv", index=False)

# 4. reported uncertainty vs empirical 1-s scatter (T1, bands 1-2)
rows = []
for d, x in data.items():
    for b in (0, 1):
        s = x["f"][:, 0, b]; u = x["u"][:, 0, b]; ok = (x["q"][:, 0, b] == 0) & (s > 0)
        m = pd.Series(np.where(ok, s, np.nan)).rolling(7, center=True, min_periods=5).median().to_numpy()
        z = ((s - m) / u)[ok]; z = z[np.isfinite(z)]
        rows.append(dict(date=d, band=b + 1, median_rel_uncertainty_pct=100 * np.median(u[ok] / s[ok]),
                         robust_std_of_scatter_over_u=1.4826 * mad(z)))
pd.DataFrame(rows).to_csv(OUT / "uncertainty_vs_scatter.csv", index=False)

# 5. scale of variability: std of log10 block means at several averaging times
rows = []
for d, x in data.items():
    for b in (0, 1):
        s = x["f"][:, 0, b]; ok = (x["q"][:, 0, b] == 0) & (s > 0)
        L = pd.Series(np.where(ok, np.log10(np.where(s > 0, s, np.nan)), np.nan))
        r = dict(date=d, band=b + 1)
        for w, n in [(1, "1s"), (60, "1min"), (600, "10min"), (1800, "30min"), (10800, "3h")]:
            r["std_log10_" + n] = L.groupby(np.arange(len(L)) // w).mean().dropna().std()
        rows.append(r)
pd.DataFrame(rows).to_csv(OUT / "variability_by_timescale.csv", index=False)

# 6. calibration drift (post-10h MAD / first-10h MAD) for three residual definitions
rows = []
for d, x in data.items():
    for b in (0, 1):
        s = pd.Series(np.where((x["q"][:, 0, b] == 0) & (x["f"][:, 0, b] > 0), x["f"][:, 0, b], np.nan))
        B = s.rolling(1800, min_periods=1800).median()
        res = {"linear": s - B, "relative": s / B - 1, "log10": np.log10(s) - np.log10(B)}
        cal = np.arange(86400) < CAL
        r = dict(date=d, band=b + 1)
        for k, v in res.items():
            v = v.to_numpy(); r[k] = mad(v[~cal & np.isfinite(v)]) / mad(v[cal & np.isfinite(v)])
        rows.append(r)
pd.DataFrame(rows).to_csv(OUT / "calibration_drift_by_residual_definition.csv", index=False)

# 7. radiation-hazard context: were any hazard channels active on these days?
rows = []
for d, x in data.items():
    ie_ok = (x["ieq"] == 0) & np.isfinite(x["ie"])
    p_ok = (x["pq"] == 0) & np.isfinite(x["p"])
    rows.append(dict(date=d, max_valid_integral_gt2MeV_electron=np.nanmax(np.where(ie_ok, x["ie"], np.nan)),
                     max_valid_proton_band4_to_11=np.nanmax(np.where(p_ok[:, :, 3:], x["p"][:, :, 3:], np.nan)),
                     max_LoLET_dose=np.nanmax(x["lo"]), max_HiLET_dose=np.nanmax(x["hi"])))
pd.DataFrame(rows).to_csv(OUT / "hazard_channel_context.csv", index=False)

pd.set_option("display.width", 200)
for f in sorted(OUT.glob("*.csv")):
    if f.name != "pct_good_by_channel.csv":
        print("\n==", f.name); print(pd.read_csv(f).round(3).to_string(index=False))
