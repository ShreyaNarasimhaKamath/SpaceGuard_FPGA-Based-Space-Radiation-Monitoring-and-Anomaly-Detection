
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy.optimize import linear_sum_assignment


# --------------------------------------------------
# CONFIGURATION
# --------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = ROOT / "data" / "raw" / "goes16_mpshi"
RESULTS_DIR = ROOT / "results" / "evaluation"
PLOTS_DIR = ROOT / "results" / "plots"

DATES = ["20240101", "20240102", "20240103"]

WINDOWS = [300, 600, 900, 1800]
REFERENCE_WINDOW = 1800

TELESCOPE_INDEX = 0
BANDS = [0, 1]  # Electron Bands 1 and 2

CALIBRATION_SAMPLES = 10 * 3600
THRESHOLD_MULTIPLIER = 6

MERGE_GAP = 120
MIN_ACTIVE_SAMPLES = 60


# --------------------------------------------------
# DATA LOADING
# --------------------------------------------------

def load_channel(date, band_index):
    filename = DATA_DIR / f"goes16_mpshi_{date}.nc"

    if not filename.exists():
        raise FileNotFoundError(f"Dataset not found: {filename}")

    with Dataset(filename, "r") as nc:
        flux_raw = nc.variables["DiffElectronFluxes"][:]
        dqf_raw = nc.variables["DiffElectronFluxDQFs"][:]

        flux = np.ma.asarray(flux_raw).astype(float).filled(np.nan)
        dqf = np.ma.asarray(dqf_raw).astype(float).filled(np.nan)

    signal = flux[:, TELESCOPE_INDEX, band_index]
    dqf = dqf[:, TELESCOPE_INDEX, band_index]

    valid = (
        np.isfinite(signal)
        & np.isfinite(dqf)
        & (dqf == 0)
        & (signal >= 0)
    )

    return signal, valid


# --------------------------------------------------
# BASELINE ESTIMATION
# --------------------------------------------------

def calculate_baseline(signal, window, method):
    series = pd.Series(signal)

    if method == "rolling_median":
        baseline = series.rolling(
            window=window,
            min_periods=window
        ).median()

    elif method == "moving_average":
        baseline = series.rolling(
            window=window,
            min_periods=window
        ).mean()

    else:
        raise ValueError(f"Unknown method: {method}")

    return baseline.to_numpy()


# --------------------------------------------------
# ROBUST THRESHOLD
# --------------------------------------------------

def calculate_threshold(residual, calibration_mask):
    calibration = residual[calibration_mask]
    calibration = calibration[np.isfinite(calibration)]

    if len(calibration) == 0:
        raise ValueError("No valid calibration residuals.")

    residual_median = np.median(calibration)
    mad = np.median(np.abs(calibration - residual_median))
    robust_sigma = 1.4826 * mad

    threshold = THRESHOLD_MULTIPLIER * robust_sigma

    return residual_median, mad, robust_sigma, threshold


# --------------------------------------------------
# EVENT GROUPING
# --------------------------------------------------

def extract_events(candidate_mask):
    indices = np.flatnonzero(candidate_mask)

    if len(indices) == 0:
        return []

    split_points = np.where(
        np.diff(indices) > MERGE_GAP
    )[0] + 1

    groups = np.split(indices, split_points)

    events = []

    for group in groups:
        active_count = len(group)

        if active_count < MIN_ACTIVE_SAMPLES:
            continue

        start = int(group[0])
        end = int(group[-1])

        events.append({
            "start": start,
            "end": end,
            "span_samples": end - start + 1,
            "active_samples": active_count,
        })

    return events


# --------------------------------------------------
# EVENT MATCHING
# --------------------------------------------------

def event_iou(event_a, event_b):
    intersection = max(
        0,
        min(event_a["end"], event_b["end"])
        - max(event_a["start"], event_b["start"]) + 1
    )

    union = (
        event_a["end"] - event_a["start"] + 1
        + event_b["end"] - event_b["start"] + 1
        - intersection
    )

    return intersection / union if union else 0.0


def match_events(reference_events, candidate_events):
    if not reference_events or not candidate_events:
        return []

    scores = np.zeros(
        (len(reference_events), len(candidate_events))
    )

    for i, ref_event in enumerate(reference_events):
        for j, candidate_event in enumerate(candidate_events):
            scores[i, j] = event_iou(ref_event, candidate_event)

    rows, cols = linear_sum_assignment(-scores)

    matches = []

    for i, j in zip(rows, cols):
        score = scores[i, j]

        if score <= 0:
            continue

        ref_event = reference_events[i]
        candidate_event = candidate_events[j]

        matches.append({
            "reference_event": i + 1,
            "candidate_event": j + 1,
            "iou": score,
            "onset_difference_seconds":
                candidate_event["start"] - ref_event["start"],
            "duration_difference_seconds":
                candidate_event["span_samples"] - ref_event["span_samples"],
        })

    return matches


# --------------------------------------------------
# MAIN EXPERIMENT
# --------------------------------------------------

def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    event_rows = []

    for date in DATES:
        print(f"\n========== DATE: {date} ==========")

        for band_index in BANDS:
            band_number = band_index + 1

            signal, valid = load_channel(date, band_index)

            if len(signal) != 86400:
                raise ValueError(
                    f"Unexpected sample count for {date}, Band {band_number}: "
                    f"{len(signal)}"
                )

            # Calculate all baselines.
            reference_baseline = calculate_baseline(
                signal, REFERENCE_WINDOW, "rolling_median"
            )

            moving_baselines = {
                window: calculate_baseline(
                    signal, window, "moving_average"
                )
                for window in WINDOWS
            }

            # Common evaluation mask:
            # all estimators must have valid baselines at each
            # timestamp used for comparison.
            common_mask = valid & np.isfinite(reference_baseline)

            for baseline in moving_baselines.values():
                common_mask &= np.isfinite(baseline)

            calibration_mask = (
                common_mask
                & (np.arange(len(signal)) < CALIBRATION_SAMPLES)
            )

            if calibration_mask.sum() < 1000:
                raise ValueError(
                    f"Insufficient common calibration samples: "
                    f"{date}, Band {band_number}"
                )

            # Reference detector.
            reference_residual = signal - reference_baseline

            _, _, _, reference_threshold = calculate_threshold(
                reference_residual, calibration_mask
            )

            reference_candidates = (
                common_mask
                & (reference_residual > reference_threshold)
            )

            reference_events = extract_events(reference_candidates)

            print(
                f"\nBand {band_number}: "
                f"valid={valid.sum()}, "
                f"common={common_mask.sum()}, "
                f"reference_events={len(reference_events)}"
            )

            for window in WINDOWS:
                baseline = moving_baselines[window]
                residual = signal - baseline

                residual_median, mad, robust_sigma, threshold = (
                    calculate_threshold(residual, calibration_mask)
                )

                candidates = (
                    common_mask
                    & (residual > threshold)
                )

                candidate_events = extract_events(candidates)

                intersection = np.logical_and(
                    reference_candidates, candidates
                ).sum()

                union = np.logical_or(
                    reference_candidates, candidates
                ).sum()

                sample_iou = intersection / union if union else 0.0

                matches = match_events(
                    reference_events, candidate_events
                )

                mean_event_iou = (
                    np.mean([match["iou"] for match in matches])
                    if matches else 0.0
                )

                summary_rows.append({
                    "date": date,
                    "band": band_number,
                    "window_samples": window,
                    "window_minutes": window / 60,
                    "common_samples": int(common_mask.sum()),
                    "calibration_samples": int(calibration_mask.sum()),
                    "threshold": threshold,
                    "residual_median": residual_median,
                    "calibration_mad": mad,
                    "robust_sigma": robust_sigma,
                    "candidate_samples": int(candidates.sum()),
                    "flagged_percentage":
                        100 * candidates.sum() / common_mask.sum(),
                    "reference_events": len(reference_events),
                    "candidate_events": len(candidate_events),
                    "matched_events": len(matches),
                    "sample_iou": sample_iou,
                    "mean_event_iou": mean_event_iou,
                })

                for match in matches:
                    event_rows.append({
                        "date": date,
                        "band": band_number,
                        "window_minutes": window / 60,
                        **match,
                    })

                print(
                    f"  MA {window // 60:2d} min | "
                    f"threshold={threshold:10.2f} | "
                    f"flagged={100 * candidates.sum() / common_mask.sum():6.2f}% | "
                    f"events={len(candidate_events):2d} | "
                    f"sample IoU={sample_iou:.3f} | "
                    f"event IoU={mean_event_iou:.3f}"
                )

    # --------------------------------------------------
    # SAVE RESULTS
    # --------------------------------------------------

    summary = pd.DataFrame(summary_rows)
    event_details = pd.DataFrame(event_rows)

    summary_path = (
        RESULTS_DIR / "moving_average_window_sensitivity.csv"
    )

    event_path = (
        RESULTS_DIR / "moving_average_window_event_matches.csv"
    )

    summary.to_csv(summary_path, index=False)
    event_details.to_csv(event_path, index=False)

    # --------------------------------------------------
    # PLOT 1: SAMPLE-LEVEL AGREEMENT
    # --------------------------------------------------

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    for band in [1, 2]:
        subset = summary[summary["band"] == band]

        means = (
            subset.groupby("window_minutes")["sample_iou"]
            .mean()
        )

        axes[0].plot(
            means.index, means.values,
            marker="o", label=f"Band {band}"
        )

    axes[0].set_title("Sample-level IoU vs Rolling Median")
    axes[0].set_xlabel("Moving-average window (minutes)")
    axes[0].set_ylabel("Mean sample IoU")
    axes[0].set_xticks([5, 10, 15, 30])
    axes[0].grid(True)
    axes[0].legend()

    # --------------------------------------------------
    # PLOT 2: FLAGGED SAMPLE PERCENTAGE
    # --------------------------------------------------

    for band in [1, 2]:
        subset = summary[summary["band"] == band]

        means = (
            subset.groupby("window_minutes")["flagged_percentage"]
            .mean()
        )

        axes[1].plot(
            means.index, means.values,
            marker="o", label=f"Band {band}"
        )

    axes[1].set_title("Mean Percentage of Samples Flagged")
    axes[1].set_xlabel("Moving-average window (minutes)")
    axes[1].set_ylabel("Flagged samples (%)")
    axes[1].set_xticks([5, 10, 15, 30])
    axes[1].grid(True)
    axes[1].legend()

    fig.tight_layout()

    plot_path = (
        PLOTS_DIR / "moving_average_window_sensitivity.png"
    )

    fig.savefig(plot_path, dpi=200)
    plt.close(fig)

    print("\nExperiment complete.")
    print(f"Summary: {summary_path}")
    print(f"Event matches: {event_path}")
    print(f"Plot: {plot_path}")


if __name__ == "__main__":
    main()
