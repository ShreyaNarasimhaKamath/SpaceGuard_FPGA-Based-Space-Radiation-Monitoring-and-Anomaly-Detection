
from pathlib import Path
import numpy as np
import pandas as pd
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parents[1]

DATA = ROOT / "data/raw/goes16_mpshi/goes16_mpshi_20240103.nc"

WINDOW = 1800
CALIBRATION_SECONDS = 36000
MULTIPLIER = 6

ALPHA = 2 / (WINDOW + 1)

EVENTS = [
    (20525, 21160),
    (37459, 38088),
    (38513, 38692),
    (57039, 57709),
    (73603, 74213),
]

METHODS = ["rolling_median", "moving_average", "ema_always"]


def mad(x):
    return np.median(np.abs(x - np.median(x)))


def main():

    with Dataset(DATA, "r") as nc:

        flux = np.ma.asarray(
            nc.variables["DiffElectronFluxes"][:, 0, 1]
        ).astype(float).filled(np.nan)

        dqf = np.ma.asarray(
            nc.variables["DiffElectronFluxDQFs"][:, 0, 1]
        ).astype(float).filled(np.nan)

    valid = (
        np.isfinite(flux)
        & np.isfinite(dqf)
        & (dqf == 0)
        & (flux >= 0)
    )

    signal = pd.Series(np.where(valid, flux, np.nan))

    baselines = {
        "rolling_median": signal.rolling(
            WINDOW, min_periods=WINDOW
        ).median(),

        "moving_average": signal.rolling(
            WINDOW, min_periods=WINDOW
        ).mean(),

        "ema_always": signal.ewm(
            alpha=ALPHA,
            adjust=False,
            min_periods=WINDOW,
            ignore_na=True
        ).mean()
    }

    results = []

    for method, baseline in baselines.items():

        residual = (signal - baseline).to_numpy()

        usable_calibration = (
            (np.arange(len(residual)) < CALIBRATION_SECONDS)
            & valid
            & np.isfinite(residual)
        )

        calibration = residual[usable_calibration]

        threshold = (
            MULTIPLIER * 1.4826 * mad(calibration)
        )

        print(f"\n{method.upper()} threshold: {threshold:.3f}")

        for event_number, (start, end) in enumerate(EVENTS, 1):

            segment = residual[start:end + 1]

            usable = np.isfinite(segment) & valid[start:end + 1]
            values = segment[usable]

            if len(values) == 0:
                continue

            active = values > threshold
            active_count = int(active.sum())

            row = {
                "event": event_number,
                "method": method,
                "start_hour": start / 3600,
                "end_hour": end / 3600,
                "max_residual": np.max(values),
                "median_residual": np.median(values),
                "threshold": threshold,
                "active_samples": active_count,
                "active_percentage": 100 * active.mean(),
                "passes_minimum_60": active_count >= 60,
            }

            results.append(row)

            print(
                f"Event {event_number}: "
                f"peak={row['max_residual']:.2f}, "
                f"median={row['median_residual']:.2f}, "
                f"active={active_count}, "
                f"active%={row['active_percentage']:.1f}"
            )

    output = pd.DataFrame(results)

    path = ROOT / "results/evaluation/jan03_band2_residuals.csv"
    output.to_csv(path, index=False)

    print(f"\nSaved: {path}")


if __name__ == "__main__":
    main()
