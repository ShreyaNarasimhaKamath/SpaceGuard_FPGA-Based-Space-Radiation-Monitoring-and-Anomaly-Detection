
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy.optimize import linear_sum_assignment


# ============================================================
# SPACEGUARD: ASYMMETRIC EMA SENSITIVITY EXPERIMENT
# ============================================================

ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = ROOT / "data" / "raw" / "goes16_mpshi"
RESULTS_DIR = ROOT / "results" / "evaluation"
PLOTS_DIR = ROOT / "results" / "plots"

DATES = ["20240101", "20240102", "20240103"]

TELESCOPE_INDEX = 0
BANDS = [0, 1]  # Electron Bands 1 and 2

REFERENCE_WINDOW = 1800
WARMUP_VALID_SAMPLES = 1800
CALIBRATION_SAMPLES = 10 * 3600

THRESHOLD_MULTIPLIER = 6
MERGE_GAP = 120
MIN_ACTIVE_SAMPLES = 60

# (alpha_up, alpha_down)
CONFIGURATIONS = {
    "A": (0.001, 0.01),
    "B": (0.002, 0.02),
    "C": (0.005, 0.02),
    "D": (0.005, 0.05),
    "E": (0.01, 0.05),
}


# ============================================================
# DATA LOADING
# ============================================================

def load_channel(date, band_index):
    filename = DATA_DIR / f"goes16_mpshi_{date}.nc"

    if not filename.exists():
        raise FileNotFoundError(filename)

    with Dataset(filename, "r") as nc:
        flux_raw = nc.variables["DiffElectronFluxes"][:]
        dqf_raw = nc.variables["DiffElectronFluxDQFs"][:]

    flux = np.ma.asarray(flux_raw).astype(float).filled(np.nan)
    dqf = np.ma.asarray(dqf_raw).astype(float).filled(np.nan)

    signal = flux[:, TELESCOPE_INDEX, band_index]
    quality = dqf[:, TELESCOPE_INDEX, band_index]

    valid = (
        np.isfinite(signal)
        & np.isfinite(quality)
        & (quality == 0)
        & (signal >= 0)
    )

    return signal, valid


# ============================================================
# ROLLING-MEDIAN REFERENCE
# ============================================================

def calculate_reference(signal):
    series = pd.Series(signal)

    return series.rolling(
        window=REFERENCE_WINDOW,
        min_periods=REFERENCE_WINDOW
    ).median().to_numpy()


# ============================================================
# ASYMMETRIC EMA
# ============================================================

def calculate_asymmetric_ema(signal, valid, alpha_up, alpha_down):
    """
    Causal asymmetric EMA.

    At each valid sample:
      1. Calculate residual against the previous baseline state.
      2. Store that baseline and residual.
      3. Update the baseline using alpha_up or alpha_down.

    Invalid samples do not update the state.

    The EMA becomes available after WARMUP_VALID_SAMPLES
    valid observations. The first valid measurement initializes
    the baseline state.
    """

    n = len(signal)

    baseline = np.full(n, np.nan)
    residual = np.full(n, np.nan)

    state = None
    valid_count = 0

    for i in range(n):
        if not valid[i]:
            continue

        x = signal[i]

        if state is None:
            state = x
            valid_count = 1
            continue

        # Use previous state to calculate the current residual.
        if valid_count >= WARMUP_VALID_SAMPLES:
            baseline[i] = state
            residual[i] = x - state

        # Update state after calculating residual.
        if x > state:
            alpha = alpha_up
        else:
            alpha = alpha_down

        state = state + alpha * (x - state)
        valid_count += 1

    return baseline, residual


# ============================================================
# ROBUST THRESHOLD
# ============================================================

def calculate_threshold(residual, calibration_mask):
    calibration = residual[calibration_mask]
    calibration = calibration[np.isfinite(calibration)]

    if len(calibration) == 0:
        raise ValueError("No valid calibration residuals.")

    median = np.median(calibration)
    mad = np.median(np.abs(calibration - median))
    robust_sigma = 1.4826 * mad

    threshold = THRESHOLD_MULTIPLIER * robust_sigma

    return median, mad, robust_sigma, threshold


# ============================================================
# EVENT EXTRACTION
# ============================================================

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


# ============================================================
# EVENT MATCHING
# ============================================================

def calculate_event_iou(event_a, event_b):
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
            scores[i, j] = calculate_event_iou(
                ref_event, candidate_event
            )

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
                candidate_event["span_samples"]
                - ref_event["span_samples"],
        })

    return matches


# ============================================================
# MAIN EXPERIMENT
# ============================================================

def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    event_rows = []

    for date in DATES:
        print(f"\n{'=' * 65}")
        print(f"DATE: {date}")
        print(f"{'=' * 65}")

        for band_index in BANDS:
            band_number = band_index + 1

            signal, valid = load_channel(date, band_index)

            if len(signal) != 86400:
                raise ValueError(
                    f"Unexpected sample count: {date}, Band {band_number}"
                )

            # Rolling-median reference.
            reference_baseline = calculate_reference(signal)
            reference_residual = signal - reference_baseline

            # Calculate all asymmetric EMA configurations.
            ema_results = {}

            for config_name, (alpha_up, alpha_down) in (
                CONFIGURATIONS.items()
            ):
                baseline, residual = calculate_asymmetric_ema(
                    signal,
                    valid,
                    alpha_up,
                    alpha_down
                )

                ema_results[config_name] = {
                    "baseline": baseline,
                    "residual": residual,
                    "alpha_up": alpha_up,
                    "alpha_down": alpha_down,
                }

            # Common evaluation timestamps across the reference
            # and all five EMA configurations.
            common_mask = (
                valid & np.isfinite(reference_residual)
            )

            for result in ema_results.values():
                common_mask &= np.isfinite(result["residual"])

            time_indices = np.arange(len(signal))

            calibration_mask = (
                common_mask
                & (time_indices < CALIBRATION_SAMPLES)
            )

            if calibration_mask.sum() < 1000:
                raise ValueError(
                    f"Insufficient calibration samples: "
                    f"{date}, Band {band_number}"
                )

            # Reference threshold and events.
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
                f"reference events={len(reference_events)}"
            )

            for config_name, result in ema_results.items():
                residual = result["residual"]

                median, mad, robust_sigma, threshold = (
                    calculate_threshold(residual, calibration_mask)
                )

                candidates = (
                    common_mask & (residual > threshold)
                )

                candidate_events = extract_events(candidates)

                intersection = np.logical_and(
                    reference_candidates, candidates
                ).sum()

                union = np.logical_or(
                    reference_candidates, candidates
                ).sum()

                sample_iou = (
                    intersection / union if union else 0.0
                )

                matches = match_events(
                    reference_events, candidate_events
                )

                mean_event_iou = (
                    np.mean([m["iou"] for m in matches])
                    if matches else 0.0
                )

                # Negative excursions are reported separately.
                negative_excursions = (
                    common_mask & (residual < -threshold)
                )

                negative_percentage = (
                    100 * negative_excursions.sum()
                    / common_mask.sum()
                )

                row = {
                    "date": date,
                    "band": band_number,
                    "configuration": config_name,
                    "alpha_up": result["alpha_up"],
                    "alpha_down": result["alpha_down"],
                    "common_samples": int(common_mask.sum()),
                    "calibration_samples": int(calibration_mask.sum()),
                    "threshold": threshold,
                    "calibration_median": median,
                    "calibration_mad": mad,
                    "robust_sigma": robust_sigma,
                    "candidate_samples": int(candidates.sum()),
                    "flagged_percentage":
                        100 * candidates.sum() / common_mask.sum(),
                    "negative_excursion_samples":
                        int(negative_excursions.sum()),
                    "negative_excursion_percentage":
                        negative_percentage,
                    "reference_events": len(reference_events),
                    "candidate_events": len(candidate_events),
                    "matched_events": len(matches),
                    "sample_iou": sample_iou,
                    "mean_event_iou": mean_event_iou,
                }

                summary_rows.append(row)

                for match in matches:
                    event_rows.append({
                        "date": date,
                        "band": band_number,
                        "configuration": config_name,
                        "alpha_up": result["alpha_up"],
                        "alpha_down": result["alpha_down"],
                        **match,
                    })

                print(
                    f"  Config {config_name} "
                    f"(up={result['alpha_up']}, "
                    f"down={result['alpha_down']}) | "
                    f"threshold={threshold:10.2f} | "
                    f"flagged={row['flagged_percentage']:6.2f}% | "
                    f"events={len(candidate_events):2d} | "
                    f"matched={len(matches):2d} | "
                    f"sample IoU={sample_iou:.3f} | "
                    f"event IoU={mean_event_iou:.3f} | "
                    f"negative={negative_percentage:.2f}%"
                )

    # ========================================================
    # SAVE CSV RESULTS
    # ========================================================

    summary = pd.DataFrame(summary_rows)
    event_details = pd.DataFrame(event_rows)

    summary_path = (
        RESULTS_DIR / "asymmetric_ema_sensitivity.csv"
    )

    event_path = (
        RESULTS_DIR / "asymmetric_ema_event_matches.csv"
    )

    summary.to_csv(summary_path, index=False)
    event_details.to_csv(event_path, index=False)

    # ========================================================
    # SUMMARY PLOTS
    # ========================================================

    fig, axes = plt.subplots(
        2, 2,
        figsize=(14, 10)
    )

    metrics = [
        ("sample_iou", "Sample-level IoU"),
        ("flagged_percentage", "Flagged samples (%)"),
        ("mean_event_iou", "Mean event-level IoU"),
        ("negative_excursion_percentage",
         "Negative excursions (%)"),
    ]

    for ax, (metric, title) in zip(axes.flat, metrics):
        for band in [1, 2]:
            subset = summary[summary["band"] == band]

            means = (
                subset.groupby("configuration")[metric]
                .mean()
                .reindex(list(CONFIGURATIONS.keys()))
            )

            ax.plot(
                means.index,
                means.values,
                marker="o",
                label=f"Band {band}"
            )

        ax.set_title(title)
        ax.set_xlabel("EMA configuration")
        ax.set_ylabel(title)
        ax.set_xticks(list(CONFIGURATIONS.keys()))
        ax.grid(True, alpha=0.3)
        ax.legend()

    fig.suptitle(
        "Asymmetric EMA Sensitivity — GOES-16 SEIS",
        fontsize=14
    )

    fig.tight_layout()

    plot_path = (
        PLOTS_DIR / "asymmetric_ema_sensitivity.png"
    )

    fig.savefig(plot_path, dpi=200)
    plt.close(fig)

    print("\n" + "=" * 65)
    print("EXPERIMENT COMPLETE")
    print("=" * 65)
    print(f"Summary CSV: {summary_path}")
    print(f"Event matches: {event_path}")
    print(f"Comparison plot: {plot_path}")


if __name__ == "__main__":
    main()
