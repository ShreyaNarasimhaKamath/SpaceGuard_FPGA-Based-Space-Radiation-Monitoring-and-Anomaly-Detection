
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = ROOT / "data/raw/goes16_mpshi"
PLOT_DIR = ROOT / "results/plots"

WINDOW = 1800
CALIBRATION_SAMPLES = 36000

# Date, Band, start hour, end hour
INTERVALS = [
    ("20240101", 1, 10, 14),
    ("20240101", 2, 10, 14),
    ("20240102", 1, 14, 18),
    ("20240102", 2, 14, 19),
    ("20240103", 1, 9, 12),
]

PLOT_DIR.mkdir(parents=True, exist_ok=True)


def main():

    for date, band, start_hour, end_hour in INTERVALS:

        filename = f"goes16_mpshi_{date}.nc"
        path = DATA_DIR / filename

        index = band - 1

        print(
            f"Inspecting {date}, Band {band}: "
            f"{start_hour}:00–{end_hour}:00"
        )

        with Dataset(path, "r") as nc:

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

        calibration_mad = np.median(
            np.abs(calibration - np.median(calibration))
        )

        threshold = 6 * 1.4826 * calibration_mad

        start = start_hour * 3600
        end = end_hour * 3600

        time = np.arange(start, end) / 3600

        fig, axes = plt.subplots(
            3, 1,
            figsize=(15, 10),
            sharex=True,
            constrained_layout=True
        )

        # Original flux and baseline.
        axes[0].plot(
            time,
            flux[start:end],
            label="Original flux",
            linewidth=1
        )

        axes[0].plot(
            time,
            baseline.iloc[start:end],
            label="30-minute rolling median",
            linewidth=1.2
        )

        axes[0].set_ylabel("Flux")
        axes[0].set_title(
            f"{date} | Band {band}: Original signal and baseline"
        )
        axes[0].legend()
        axes[0].grid(alpha=0.3)

        # Residual.
        r = residual.iloc[start:end].to_numpy()

        axes[1].plot(
            time,
            r,
            linewidth=0.8,
            label="Residual"
        )

        axes[1].axhline(
            threshold,
            linestyle="--",
            label="+6 sigma"
        )

        axes[1].axhline(
            -threshold,
            linestyle="--",
            label="-6 sigma"
        )

        axes[1].axhline(0, color="black", linewidth=0.7)

        axes[1].set_ylabel("Residual flux")
        axes[1].legend()
        axes[1].grid(alpha=0.3)

        # Validity.
        axes[2].step(
            time,
            valid[start:end].astype(int),
            where="post"
        )

        axes[2].set_ylabel("Input valid")
        axes[2].set_xlabel("Hour of day")
        axes[2].set_ylim(-0.1, 1.1)
        axes[2].grid(alpha=0.3)

        fig.savefig(
            PLOT_DIR /
            f"baseline_transition_{date}_band{band}.png",
            dpi=180
        )

        plt.close(fig)

    print("\nBASELINE TRANSITION INSPECTION COMPLETE")


if __name__ == "__main__":
    main()
