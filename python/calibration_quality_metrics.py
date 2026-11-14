
from pathlib import Path
import numpy as np
import pandas as pd
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/raw/goes16_mpshi"

FILES = [
    "goes16_mpshi_20240101.nc",
    "goes16_mpshi_20240102.nc",
    "goes16_mpshi_20240103.nc",
]

WINDOW = 1800
CALIBRATION_SAMPLES = 36000


def mad(values):
    center = np.median(values)
    return np.median(np.abs(values - center))


rows = []

for filename in FILES:
    date = filename.split("_")[-1].split(".")[0]
    path = DATA / filename

    with Dataset(path, "r") as nc:
        for band, index in [(1, 0), (2, 1)]:

            flux = np.ma.asarray(
                nc.variables["DiffElectronFluxes"][:, 0, index]
            ).astype(float).filled(np.nan)

            dqf = np.ma.asarray(
                nc.variables["DiffElectronFluxDQFs"][:, 0, index]
            ).astype(np.int16).filled(-1)

            valid = (
                (dqf == 0)
                & np.isfinite(flux)
                & (flux >= 0)
            )

            signal = pd.Series(
                np.where(valid, flux, np.nan)
            )

            baseline = signal.rolling(
                WINDOW,
                min_periods=WINDOW
            ).median()

            residual = signal - baseline

            calibration = residual.iloc[
                :CALIBRATION_SAMPLES
            ].dropna().to_numpy()

            post = residual.iloc[
                CALIBRATION_SAMPLES:
            ].dropna().to_numpy()

            calibration_mad = mad(calibration)
            post_mad = mad(post)

            calibration_sigma = 1.4826 * calibration_mad
            post_sigma = 1.4826 * post_mad

            rows.append({
                "date": date,
                "band": band,

                "input_valid_samples": int(valid.sum()),
                "input_invalid_samples": int((~valid).sum()),

                "baseline_available_samples": int(
                    baseline.notna().sum()
                ),

                "baseline_unavailable_samples": int(
                    baseline.isna().sum()
                ),

                "calibration_residual_samples": len(calibration),
                "post_residual_samples": len(post),

                "calibration_median": float(np.median(calibration)),
                "post_median": float(np.median(post)),

                "calibration_mad": float(calibration_mad),
                "post_mad": float(post_mad),

                "calibration_robust_sigma": float(
                    calibration_sigma
                ),

                "post_robust_sigma": float(post_sigma),

                "post_to_calibration_mad_ratio": float(
                    post_mad / calibration_mad
                    if calibration_mad > 0 else np.nan
                ),

                "post_to_calibration_sigma_ratio": float(
                    post_sigma / calibration_sigma
                    if calibration_sigma > 0 else np.nan
                ),
            })

df = pd.DataFrame(rows)

output = ROOT / "results/evaluation/calibration_quality_metrics.csv"
output.parent.mkdir(parents=True, exist_ok=True)

df.to_csv(output, index=False)

print("\nCALIBRATION QUALITY METRICS")
print(df.to_string(index=False))

print("\nSaved:", output)
