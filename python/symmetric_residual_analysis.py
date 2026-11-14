
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = ROOT / "data/raw/goes16_mpshi"
PLOT_DIR = ROOT / "results/plots"
EVAL_DIR = ROOT / "results/evaluation"

FILES = [
    "goes16_mpshi_20240101.nc",
    "goes16_mpshi_20240102.nc",
    "goes16_mpshi_20240103.nc",
]

WINDOW = 1800
CALIBRATION_SAMPLES = 36000
MULTIPLIER = 6

PLOT_DIR.mkdir(parents=True, exist_ok=True)
EVAL_DIR.mkdir(parents=True, exist_ok=True)


def mad(values):
    if len(values) == 0:
        return np.nan

    center = np.median(values)
    return np.median(np.abs(values - center))


records = []

for filename in FILES:

    date = filename.split("_")[-1].split(".")[0]
    print(f"\nProcessing {date}")

    with Dataset(DATA_DIR / filename, "r") as nc:

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

            residual = (signal - baseline).to_numpy()

            calibration = residual[:CALIBRATION_SAMPLES]
            calibration = calibration[np.isfinite(calibration)]

            calibration_mad = mad(calibration)
            threshold = MULTIPLIER * 1.4826 * calibration_mad

            for hour in range(24):

                start = hour * 3600
                end = (hour + 1) * 3600

                values = residual[start:end]
                values = values[np.isfinite(values)]

                if len(values) == 0:
                    continue

                positive = values > threshold
                negative = values < -threshold
                central = ~(positive | negative)

                central_values = values[central]

                records.append({

                    "date": date,
                    "band": band,
                    "hour": hour,

                    "valid_samples": len(values),

                    "all_mad": mad(values),

                    "central_samples": len(central_values),
                    "central_mad": mad(central_values),

                    "positive_excursion_percent":
                        100 * positive.mean(),

                    "negative_excursion_percent":
                        100 * negative.mean(),

                    "median": np.median(values),

                    "p05": np.percentile(values, 5),
                    "p95": np.percentile(values, 95),

                    "threshold": threshold,

                    "all_mad_ratio":
                        mad(values) / calibration_mad,

                    "central_mad_ratio":
                        (
                            mad(central_values) / calibration_mad
                            if len(central_values) > 0
                            else np.nan
                        )
                })

df = pd.DataFrame(records)

csv_path = EVAL_DIR / "symmetric_residual_analysis.csv"
df.to_csv(csv_path, index=False)

# Plot all MAD against central MAD.
fig, axes = plt.subplots(
    3, 2, figsize=(15, 12), constrained_layout=True
)

for row, filename in enumerate(FILES):

    date = filename.split("_")[-1].split(".")[0]

    for col, band in enumerate([1, 2]):

        ax = axes[row, col]

        subset = df[
            (df["date"] == date)
            & (df["band"] == band)
        ]

        ax.plot(
            subset["hour"],
            subset["all_mad_ratio"],
            marker="o",
            label="All residuals"
        )

        ax.plot(
            subset["hour"],
            subset["central_mad_ratio"],
            marker="s",
            label="Central residuals"
        )

        ax.axhline(1, linestyle="--")
        ax.axvline(10, linestyle=":")

        ax.set_title(f"{date} | Band {band}")
        ax.set_xlabel("Hour")
        ax.set_ylabel("MAD / calibration MAD")
        ax.set_xticks(range(0, 24, 2))
        ax.grid(True, alpha=0.3)
        ax.legend()

fig.suptitle("Symmetric Residual Scale Analysis")

fig.savefig(
    PLOT_DIR / "symmetric_residual_scale.png",
    dpi=200
)

plt.close(fig)

# Plot positive and negative excursion percentages.
fig, axes = plt.subplots(
    3, 2, figsize=(15, 12), constrained_layout=True
)

for row, filename in enumerate(FILES):

    date = filename.split("_")[-1].split(".")[0]

    for col, band in enumerate([1, 2]):

        ax = axes[row, col]

        subset = df[
            (df["date"] == date)
            & (df["band"] == band)
        ]

        ax.plot(
            subset["hour"],
            subset["positive_excursion_percent"],
            marker="o",
            label="Positive excursions"
        )

        ax.plot(
            subset["hour"],
            subset["negative_excursion_percent"],
            marker="s",
            label="Negative excursions"
        )

        ax.axvline(10, linestyle=":")

        ax.set_title(f"{date} | Band {band}")
        ax.set_xlabel("Hour")
        ax.set_ylabel("Excursion samples (%)")
        ax.set_xticks(range(0, 24, 2))
        ax.grid(True, alpha=0.3)
        ax.legend()

fig.suptitle("Positive vs Negative Residual Excursions")

fig.savefig(
    PLOT_DIR / "positive_negative_excursions.png",
    dpi=200
)

plt.close(fig)

print("\nSYMMETRIC RESIDUAL ANALYSIS COMPLETE")
print("\nSaved:", csv_path)
print(PLOT_DIR / "symmetric_residual_scale.png")
print(PLOT_DIR / "positive_negative_excursions.png")

print("\nCentral MAD ratio summary:")

print(
    df.groupby(["date", "band"])[
        "central_mad_ratio"
    ].agg(["median", "max"]).to_string()
)
