
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

BANDS = {1: 0, 2: 1}

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


def main():

    records = []

    for filename in FILES:

        date = filename.split("_")[-1].split(".")[0]
        path = DATA_DIR / filename

        print(f"\nProcessing {date}")

        with Dataset(path, "r") as nc:

            for band, index in BANDS.items():

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
                threshold = (
                    MULTIPLIER * 1.4826 * calibration_mad
                )

                for hour in range(24):

                    start = hour * 3600
                    end = (hour + 1) * 3600

                    hour_residual = residual[start:end]
                    hour_valid = valid[start:end]

                    usable = (
                        np.isfinite(hour_residual)
                        & hour_valid
                    )

                    values = hour_residual[usable]

                    if len(values) == 0:
                        continue

                    candidates = values > threshold
                    non_candidates = values[~candidates]

                    records.append({

                        "date": date,
                        "band": band,
                        "hour": hour,

                        "valid_residual_samples": len(values),

                        "candidate_samples": int(candidates.sum()),

                        "candidate_percent": (
                            100 * candidates.mean()
                        ),

                        "all_residual_median": float(
                            np.median(values)
                        ),

                        "all_residual_mad": float(
                            mad(values)
                        ),

                        "below_threshold_samples": len(
                            non_candidates
                        ),

                        "below_threshold_median": (
                            float(np.median(non_candidates))
                            if len(non_candidates) else np.nan
                        ),

                        "below_threshold_mad": mad(
                            non_candidates
                        ),

                        "calibration_mad": float(
                            calibration_mad
                        ),

                        "all_mad_ratio": (
                            mad(values) / calibration_mad
                            if calibration_mad > 0 else np.nan
                        ),

                        "below_threshold_mad_ratio": (
                            mad(non_candidates) / calibration_mad
                            if len(non_candidates) and calibration_mad > 0
                            else np.nan
                        ),
                    })

    df = pd.DataFrame(records)

    csv_path = (
        EVAL_DIR / "background_candidate_separation.csv"
    )

    df.to_csv(csv_path, index=False)

    # --------------------------------------------------------
    # Plot MAD comparison
    # --------------------------------------------------------

    fig, axes = plt.subplots(
        3, 2,
        figsize=(15, 12),
        constrained_layout=True
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
                label="All residual samples"
            )

            ax.plot(
                subset["hour"],
                subset["below_threshold_mad_ratio"],
                marker="s",
                label="Below 6x threshold"
            )

            ax.axhline(1, linestyle="--")
            ax.axvline(10, linestyle=":")

            ax.set_title(f"{date} | Band {band}")
            ax.set_xlabel("Hour")
            ax.set_ylabel("MAD / calibration MAD")
            ax.set_xticks(range(0, 24, 2))
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=8)

    fig.suptitle(
        "Residual Scale: All Samples vs Below-Threshold Samples"
    )

    fig.savefig(
        PLOT_DIR / "background_candidate_separation.png",
        dpi=200
    )

    plt.close(fig)

    # --------------------------------------------------------
    # Plot candidate percentages
    # --------------------------------------------------------

    fig, axes = plt.subplots(
        3, 2,
        figsize=(15, 12),
        constrained_layout=True
    )

    for row, filename in enumerate(FILES):

        date = filename.split("_")[-1].split(".")[0]

        for col, band in enumerate([1, 2]):

            ax = axes[row, col]

            subset = df[
                (df["date"] == date)
                & (df["band"] == band)
            ]

            ax.bar(
                subset["hour"],
                subset["candidate_percent"]
            )

            ax.axvline(10, linestyle=":")
            ax.set_title(f"{date} | Band {band}")
            ax.set_xlabel("Hour")
            ax.set_ylabel("Candidate samples (%)")
            ax.set_xticks(range(0, 24, 2))
            ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle("Hourly Candidate Activity")

    fig.savefig(
        PLOT_DIR / "candidate_activity_by_hour.png",
        dpi=200
    )

    plt.close(fig)

    print("\nANALYSIS COMPLETE")
    print("\nSaved:")
    print(csv_path)
    print(PLOT_DIR / "background_candidate_separation.png")
    print(PLOT_DIR / "candidate_activity_by_hour.png")


if __name__ == "__main__":
    main()
