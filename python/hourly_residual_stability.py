
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
    center = np.median(values)
    return np.median(np.abs(values - center))


def main():

    records = []
    plot_data = {}

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

                residual = signal - baseline

                calibration = residual.iloc[
                    :CALIBRATION_SAMPLES
                ].dropna().to_numpy()

                center = np.median(calibration)
                calibration_mad = mad(calibration)

                robust_sigma = 1.4826 * calibration_mad
                threshold = MULTIPLIER * robust_sigma

                hourly_mad = []
                hourly_sigma = []

                for hour in range(24):

                    start = hour * 3600
                    end = (hour + 1) * 3600

                    hour_residual = residual.iloc[
                        start:end
                    ].dropna().to_numpy()

                    hour_valid = valid[start:end]

                    hour_invalid = int(
                        np.sum(~hour_valid)
                    )

                    if len(hour_residual) == 0:
                        hour_mad_value = np.nan
                        hour_sigma_value = np.nan
                        median_value = np.nan
                        flagged_percent = np.nan

                    else:
                        hour_mad_value = mad(hour_residual)

                        hour_sigma_value = (
                            1.4826 * hour_mad_value
                        )

                        median_value = float(
                            np.median(hour_residual)
                        )

                        flagged_percent = (
                            100.0
                            * np.sum(hour_residual > threshold)
                            / len(hour_residual)
                        )

                    hourly_mad.append(hour_mad_value)
                    hourly_sigma.append(hour_sigma_value)

                    records.append({

                        "date": date,
                        "band": band,
                        "hour": hour,

                        "valid_residual_samples": len(
                            hour_residual
                        ),

                        "invalid_input_samples": hour_invalid,

                        "residual_median": median_value,
                        "hourly_mad": hour_mad_value,
                        "hourly_robust_sigma": hour_sigma_value,

                        "calibration_mad": float(
                            calibration_mad
                        ),

                        "calibration_robust_sigma": float(
                            robust_sigma
                        ),

                        "threshold_6x": float(threshold),

                        "hourly_to_calibration_mad_ratio": (
                            hour_mad_value / calibration_mad
                            if calibration_mad > 0 else np.nan
                        ),

                        "flagged_percent": flagged_percent,
                    })

                plot_data[(date, band)] = {
                    "mad": hourly_mad,
                    "sigma": hourly_sigma,
                    "threshold": threshold,
                    "calibration_mad": calibration_mad,
                }

                print(
                    f"Band {band}: calibration MAD="
                    f"{calibration_mad:.3f}, "
                    f"6x threshold={threshold:.3f}"
                )

    # Save hourly results.
    df = pd.DataFrame(records)

    output_csv = (
        EVAL_DIR / "hourly_residual_stability.csv"
    )

    df.to_csv(output_csv, index=False)

    # --------------------------------------------------------
    # Plot hourly MAD ratios
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

            data = plot_data[(date, band)]

            ratios = (
                np.array(data["mad"])
                / data["calibration_mad"]
            )

            ax.plot(
                range(24),
                ratios,
                marker="o"
            )

            ax.axhline(
                1.0,
                linestyle="--",
                label="Calibration MAD"
            )

            ax.axvline(
                10,
                linestyle=":",
                label="Calibration boundary"
            )

            ax.set_title(f"{date} | Band {band}")
            ax.set_xlabel("Hour")
            ax.set_ylabel("Hourly MAD / calibration MAD")
            ax.set_xticks(range(0, 24, 2))
            ax.grid(True, alpha=0.3)
            ax.legend()

    fig.suptitle("Hourly Residual Scale Stability")

    fig.savefig(
        PLOT_DIR / "hourly_residual_stability.png",
        dpi=200
    )

    plt.close(fig)

    # --------------------------------------------------------
    # Plot hourly flagged percentages
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
                subset["flagged_percent"],
                marker="o"
            )

            ax.axvline(
                10,
                linestyle=":",
                label="Calibration boundary"
            )

            ax.set_title(f"{date} | Band {band}")
            ax.set_xlabel("Hour")
            ax.set_ylabel("Samples above 6x threshold (%)")
            ax.set_xticks(range(0, 24, 2))
            ax.grid(True, alpha=0.3)
            ax.legend()

    fig.suptitle("Hourly Candidate Sample Percentage")

    fig.savefig(
        PLOT_DIR / "hourly_flagged_percentage.png",
        dpi=200
    )

    plt.close(fig)

    print("\nHOURLY ANALYSIS COMPLETE")
    print("\nSaved:")
    print(output_csv)
    print(PLOT_DIR / "hourly_residual_stability.png")
    print(PLOT_DIR / "hourly_flagged_percentage.png")

    print("\nHourly MAD ratio summary:")
    print(
        df.groupby(["date", "band"])[
            "hourly_to_calibration_mad_ratio"
        ].agg(["median", "max"]).to_string()
    )


if __name__ == "__main__":
    main()
