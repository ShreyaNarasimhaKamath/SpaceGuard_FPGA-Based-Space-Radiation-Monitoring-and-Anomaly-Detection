
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
MULTIPLIERS = [4, 5, 6, 7, 8]

WINDOW = 1800
CALIBRATION_SAMPLES = 10 * 3600
REFERENCE_MULTIPLIER = 6
MERGE_GAP = 120
MIN_ACTIVE = 60

PLOT_DIR.mkdir(parents=True, exist_ok=True)
EVAL_DIR.mkdir(parents=True, exist_ok=True)


def load_channel(path, band_index):
    with Dataset(path, "r") as nc:
        flux = np.ma.asarray(
            nc.variables["DiffElectronFluxes"][
                :, 0, band_index
            ]
        ).astype(float).filled(np.nan)

        dqf = np.ma.asarray(
            nc.variables["DiffElectronFluxDQFs"][
                :, 0, band_index
            ]
        ).astype(np.int16).filled(-1)

    valid = (
        (dqf == 0)
        & np.isfinite(flux)
        & (flux >= 0)
    )

    signal = pd.Series(
        np.where(valid, flux, np.nan)
    )

    return flux, valid, signal


def detect(signal, valid, flux, multiplier):
    baseline = signal.rolling(
        WINDOW,
        min_periods=WINDOW
    ).median()

    residual = signal - baseline

    calibration = residual.iloc[
        :CALIBRATION_SAMPLES
    ].dropna()

    if calibration.empty:
        raise ValueError("Calibration residual is empty")

    center = calibration.median()

    mad = np.median(
        np.abs(calibration.to_numpy() - center)
    )

    sigma = 1.4826 * mad
    threshold = multiplier * sigma

    candidate = (
        (residual > threshold)
        & valid
    ).to_numpy()

    indices = np.flatnonzero(candidate)
    events = []

    if len(indices):
        start = previous = indices[0]
        active = 1

        groups = []

        for index in indices[1:]:
            if index - previous - 1 <= MERGE_GAP:
                active += 1
            else:
                groups.append((start, previous, active))
                start = index
                active = 1

            previous = index

        groups.append((start, previous, active))

        for start, end, active in groups:
            if active < MIN_ACTIVE:
                continue

            active_indices = indices[
                (indices >= start) & (indices <= end)
            ]

            events.append({
                "start_index": int(start),
                "end_index": int(end),
                "start_hour": start / 3600,
                "end_hour": end / 3600,
                "active_samples": int(active),
                "peak_flux": float(np.max(flux[active_indices])),
            })

    return residual, threshold, candidate, events, sigma


def main():
    residual_records = []
    event_records = []
    residual_data = {}

    for filename in FILES:
        date = filename.split("_")[-1].split(".")[0].replace("2024", "2024")

        print(f"\nProcessing {date}")

        path = DATA_DIR / filename

        for band, band_index in BANDS.items():
            flux, valid, signal = load_channel(path, band_index)

            # Use the existing reference threshold for distributions.
            residual, threshold, candidate, events, sigma = detect(
                signal, valid, flux, REFERENCE_MULTIPLIER
            )

            calibration = residual.iloc[:CALIBRATION_SAMPLES].dropna()
            valid_residual = residual.dropna()

            residual_data[(date, band)] = {
                "residual": valid_residual.to_numpy(),
                "threshold": threshold,
            }

            residual_records.append({
                "date": date,
                "band": band,
                "valid_samples": int(valid.sum()),
                "calibration_residual_samples": len(calibration),
                "residual_median": float(valid_residual.median()),
                "residual_mad_calibration": float(
                    np.median(np.abs(calibration.to_numpy() - calibration.median()))
                ),
                "robust_sigma": float(sigma),
                "threshold_6x": float(threshold),
                "residual_std": float(valid_residual.std()),
                "residual_95th_percentile": float(valid_residual.quantile(0.95)),
                "residual_99th_percentile": float(valid_residual.quantile(0.99)),
                "candidate_samples": int(candidate.sum()),
                "event_count": len(events),
            })

            for event_id, event in enumerate(events, start=1):
                event_records.append({
                    "date": date,
                    "band": band,
                    "event_id": event_id,
                    **event,
                })

    # --------------------------------------------------------
    # Residual distribution plots
    # --------------------------------------------------------

    fig, axes = plt.subplots(
        3, 2, figsize=(14, 13), constrained_layout=True
    )

    for row, filename in enumerate(FILES):
        date = filename.split("_")[-1].split(".")[0]

        for col, band in enumerate([1, 2]):
            ax = axes[row, col]
            data = residual_data[(date, band)]

            values = data["residual"]
            threshold = data["threshold"]

            ax.hist(
                values,
                bins=120,
                density=True,
                alpha=0.75
            )

            ax.axvline(
                threshold,
                linestyle="--",
                label=f"+6σ threshold = {threshold:.1f}"
            )

            ax.set_title(f"{date} | Band {band}")
            ax.set_xlabel("Residual flux")
            ax.set_ylabel("Density")
            ax.legend()

    fig.suptitle("GOES-16 Residual Distributions (6× threshold)")
    fig.savefig(
        PLOT_DIR / "residual_distributions.png",
        dpi=200
    )
    plt.close(fig)

    # --------------------------------------------------------
    # Threshold comparison plot
    # --------------------------------------------------------

    summary = pd.DataFrame(residual_records)

    fig, axes = plt.subplots(
        3, 2, figsize=(13, 11), constrained_layout=True
    )

    for row, filename in enumerate(FILES):
        date = filename.split("_")[-1].split(".")[0]

        for col, band in enumerate([1, 2]):
            ax = axes[row, col]

            record = summary[
                (summary["date"] == date)
                & (summary["band"] == band)
            ].iloc[0]

            # Recompute detection counts at each multiplier.
            path = DATA_DIR / filename
            flux, valid, signal = load_channel(
                path, BANDS[band]
            )

            flagged = []

            for multiplier in MULTIPLIERS:
                _, _, candidate, _, _ = detect(
                    signal, valid, flux, multiplier
                )

                flagged.append(
                    100 * candidate.sum() / valid.sum()
                )

            ax.plot(
                MULTIPLIERS,
                flagged,
                marker="o"
            )

            ax.set_title(f"{date} | Band {band}")
            ax.set_xlabel("Threshold multiplier")
            ax.set_ylabel("Flagged samples (%)")
            ax.set_xticks(MULTIPLIERS)
            ax.grid(True, alpha=0.3)

    fig.suptitle("Threshold Sensitivity by Date and Energy Band")
    fig.savefig(
        PLOT_DIR / "threshold_comparison.png",
        dpi=200
    )
    plt.close(fig)

    # --------------------------------------------------------
    # Event persistence across threshold settings
    # --------------------------------------------------------

    for filename in FILES:
        date = filename.split("_")[-1].split(".")[0]
        path = DATA_DIR / filename

        for band, band_index in BANDS.items():
            flux, valid, signal = load_channel(path, band_index)

            detections = {}

            for multiplier in MULTIPLIERS:
                _, _, _, events, _ = detect(
                    signal, valid, flux, multiplier
                )
                detections[multiplier] = events

            reference_events = detections[REFERENCE_MULTIPLIER]

            for event_id, reference in enumerate(
                reference_events, start=1
            ):
                row = {
                    "date": date,
                    "band": band,
                    "reference_multiplier": REFERENCE_MULTIPLIER,
                    "reference_event_id": event_id,
                    "reference_start_hour": reference["start_hour"],
                    "reference_end_hour": reference["end_hour"],
                }

                for multiplier in MULTIPLIERS:
                    best_iou = 0.0
                    best_event = None

                    for candidate_event in detections[multiplier]:
                        intersection = max(
                            0,
                            min(
                                reference["end_index"],
                                candidate_event["end_index"]
                            )
                            - max(
                                reference["start_index"],
                                candidate_event["start_index"]
                            ) + 1
                        )

                        union = (
                            reference["end_index"]
                            - reference["start_index"] + 1
                            + candidate_event["end_index"]
                            - candidate_event["start_index"] + 1
                            - intersection
                        )

                        iou = intersection / union if union else 0

                        if iou > best_iou:
                            best_iou = iou
                            best_event = candidate_event

                    row[f"multiplier_{multiplier}_best_iou"] = best_iou
                    row[f"multiplier_{multiplier}_matched"] = (
                        best_event is not None
                    )

                event_records.append(row)

    # Save residual statistics.
    summary.to_csv(
        EVAL_DIR / "residual_statistics.csv",
        index=False
    )

    # Event persistence has its own output file.
    persistence_df = pd.DataFrame(event_records)

    persistence_df.to_csv(
        EVAL_DIR / "event_persistence.csv",
        index=False
    )

    print("\nANALYSIS COMPLETE")
    print("\nResidual statistics:")
    print(summary.to_string(index=False))

    print("\nGenerated files:")
    print(PLOT_DIR / "residual_distributions.png")
    print(PLOT_DIR / "threshold_comparison.png")
    print(EVAL_DIR / "residual_statistics.csv")
    print(EVAL_DIR / "event_persistence.csv")


if __name__ == "__main__":
    main()

