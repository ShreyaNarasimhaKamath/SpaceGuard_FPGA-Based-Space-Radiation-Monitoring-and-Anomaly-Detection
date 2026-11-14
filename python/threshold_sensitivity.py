
from pathlib import Path
import re

import numpy as np
import pandas as pd
from netCDF4 import Dataset
from scipy.optimize import linear_sum_assignment


# ============================================================
# CONFIGURATION — KEEP THESE FIXED FOR THE EXPERIMENT
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = PROJECT_ROOT / "data" / "raw" / "goes16_mpshi"
OUTPUT_DIR = PROJECT_ROOT / "results" / "evaluation"

FILES = [
    "goes16_mpshi_20240101.nc",
    "goes16_mpshi_20240102.nc",
    "goes16_mpshi_20240103.nc",
]

THRESHOLD_MULTIPLIERS = [4, 5, 6, 7, 8]

TELESCOPE_INDEX = 0
BANDS = {
    1: 0,
    2: 1,
}

BASELINE_WINDOW = 1800
CALIBRATION_HOURS = 10
MERGE_GAP_SECONDS = 120
MIN_ACTIVE_SAMPLES = 60

FLUX_VARIABLE = "DiffElectronFluxes"
DQF_VARIABLE = "DiffElectronFluxDQFs"


# ============================================================
# DATA LOADING
# ============================================================


def load_channel(file_path, band_index):
    """
    Load one telescope/energy-band channel.

    A sample is valid only when:
      - DQF == 0
      - flux is finite
      - flux >= 0
    """

    with Dataset(file_path, "r") as nc:
        flux_data = nc.variables[FLUX_VARIABLE][
            :, TELESCOPE_INDEX, band_index
        ]

        dqf_data = nc.variables[DQF_VARIABLE][
            :, TELESCOPE_INDEX, band_index
        ]

        # Convert DQF to signed integer before filling masked values.
        flux = np.ma.asarray(flux_data).astype(float).filled(np.nan)

        dqf = (
            np.ma.asarray(dqf_data)
            .astype(np.int16)
            .filled(-1)
        )

    valid = (
        (dqf == 0)
        & np.isfinite(flux)
        & (flux >= 0)
    )

    signal = pd.Series(
        np.where(valid, flux, np.nan)
    )

    return flux, valid, signal


# ============================================================
# DETECTION
# ============================================================

def detect_events(flux, valid, signal, multiplier):
    """
    Apply the causal rolling-median detector.

    Threshold is calibrated independently for each channel
    using the first 10 hours of residual samples.
    """

    baseline = signal.rolling(
        window=BASELINE_WINDOW,
        min_periods=BASELINE_WINDOW
    ).median()

    residual = signal - baseline

    calibration_samples = int(CALIBRATION_HOURS * 3600)

    calibration = residual.iloc[:calibration_samples].dropna()

    if len(calibration) == 0:
        raise ValueError("No valid calibration residual samples.")

    residual_median = calibration.median()

    mad = np.median(
        np.abs(calibration.to_numpy() - residual_median)
    )

    robust_sigma = 1.4826 * mad
    threshold = multiplier * robust_sigma

    if not np.isfinite(threshold) or threshold <= 0:
        raise ValueError(
            f"Invalid threshold: {threshold}. "
            "Check calibration data."
        )

    candidate = (
        (residual > threshold)
        & valid
    ).to_numpy()

    candidate_indices = np.flatnonzero(candidate)

    events = []

    if len(candidate_indices) == 0:
        return {
            "threshold": threshold,
            "robust_sigma": robust_sigma,
            "candidate_count": 0,
            "candidate_percent": 0.0,
            "events": events,
        }

    # Group candidate samples into events.
    # A gap of at most MERGE_GAP_SECONDS is merged.
    groups = []

    group_start = candidate_indices[0]
    previous_index = candidate_indices[0]
    active_count = 1

    for index in candidate_indices[1:]:
        gap = index - previous_index - 1

        if gap <= MERGE_GAP_SECONDS:
            active_count += 1
        else:
            groups.append(
                (group_start, previous_index, active_count)
            )

            group_start = index
            active_count = 1

        previous_index = index

    groups.append(
        (group_start, previous_index, active_count)
    )

    # Retain events meeting the minimum active-sample requirement.
    for start, end, active_count in groups:

        if active_count < MIN_ACTIVE_SAMPLES:
            continue

        active_in_event = candidate[start:end + 1]
        active_indices = np.flatnonzero(active_in_event) + start

        peak_flux = float(np.max(flux[active_indices]))

        events.append({
            "start_index": int(start),
            "end_index": int(end),
            "start_hour": start / 3600.0,
            "end_hour": end / 3600.0,
            "span_duration_seconds": int(end - start + 1),
            "active_samples": int(active_count),
            "peak_flux": peak_flux,
        })

    valid_count = int(np.sum(valid))
    candidate_count = int(np.sum(candidate))

    candidate_percent = (
        100.0 * candidate_count / valid_count
        if valid_count else 0.0
    )

    return {
        "threshold": float(threshold),
        "robust_sigma": float(robust_sigma),
        "candidate_count": candidate_count,
        "candidate_percent": candidate_percent,
        "events": events,
    }


# ============================================================
# CROSS-CHANNEL EVENT MATCHING
# ============================================================

def calculate_iou(event_a, event_b):
    """
    Calculate temporal intersection over union.

    Event intervals are treated as inclusive sample intervals.
    """

    start_a = event_a["start_index"]
    end_a = event_a["end_index"]

    start_b = event_b["start_index"]
    end_b = event_b["end_index"]

    intersection = max(
        0,
        min(end_a, end_b) - max(start_a, start_b) + 1
    )

    if intersection == 0:
        return 0.0

    duration_a = end_a - start_a + 1
    duration_b = end_b - start_b + 1

    union = duration_a + duration_b - intersection

    return intersection / union


def match_events(band1_events, band2_events):
    """
    One-to-one maximum-IoU assignment.
    Only positive-overlap pairs are retained.
    """

    if not band1_events or not band2_events:
        return [], 0.0

    iou_matrix = np.zeros(
        (len(band1_events), len(band2_events))
    )

    for i, event1 in enumerate(band1_events):
        for j, event2 in enumerate(band2_events):
            iou_matrix[i, j] = calculate_iou(event1, event2)

    row_indices, col_indices = linear_sum_assignment(
        -iou_matrix
    )

    matches = []

    for i, j in zip(row_indices, col_indices):
        iou = float(iou_matrix[i, j])

        if iou <= 0:
            continue

        event1 = band1_events[i]
        event2 = band2_events[j]

        overlap_seconds = max(
            0,
            min(
                event1["end_index"],
                event2["end_index"]
            )
            - max(
                event1["start_index"],
                event2["start_index"]
            )
            + 1
        )

        matches.append({
            "band1_event": i + 1,
            "band2_event": j + 1,
            "iou": iou,
            "overlap_seconds": int(overlap_seconds),
        })

    mean_iou = (
        float(np.mean([m["iou"] for m in matches]))
        if matches else 0.0
    )

    return matches, mean_iou


# ============================================================
# MAIN EXPERIMENT
# ============================================================

def main():

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    event_rows = []
    match_rows = []

    for filename in FILES:

        file_path = DATA_DIR / filename

        if not file_path.exists():
            raise FileNotFoundError(
                f"Dataset file not found: {file_path}"
            )

        date_match = re.search(r"(\d{8})", filename)

        if not date_match:
            raise ValueError(
                f"Could not extract date from filename: {filename}"
            )

        date = date_match.group(1)

        print("\n" + "=" * 65)
        print(f"PROCESSING DATE: {date}")
        print("=" * 65)

        # Load both channels once per day.
        channel_data = {}

        for band_number, band_index in BANDS.items():

            flux, valid, signal = load_channel(
                file_path,
                band_index
            )

            channel_data[band_number] = {
                "flux": flux,
                "valid": valid,
                "signal": signal,
            }

        # Repeat detection for each threshold multiplier.
        for multiplier in THRESHOLD_MULTIPLIERS:

            print(f"\n--- Threshold multiplier: {multiplier} ---")

            detected = {}

            for band_number in BANDS:

                data = channel_data[band_number]

                result = detect_events(
                    data["flux"],
                    data["valid"],
                    data["signal"],
                    multiplier
                )

                detected[band_number] = result

                print(
                    f"Band {band_number}: "
                    f"threshold={result['threshold']:.3f}, "
                    f"candidates={result['candidate_count']}, "
                    f"flagged={result['candidate_percent']:.3f}%, "
                    f"events={len(result['events'])}"
                )

                # Save individual event records.
                for event_number, event in enumerate(
                    result["events"], start=1
                ):

                    event_rows.append({
                        "date": date,
                        "threshold_multiplier": multiplier,
                        "band": band_number,
                        "event_id": event_number,
                        **event,
                    })

            # Match Band 1 and Band 2 events.
            band1_events = detected[1]["events"]
            band2_events = detected[2]["events"]

            matches, mean_iou = match_events(
                band1_events,
                band2_events
            )

            matched_b1 = {m["band1_event"] for m in matches}
            matched_b2 = {m["band2_event"] for m in matches}

            unmatched_b1 = len(band1_events) - len(matched_b1)
            unmatched_b2 = len(band2_events) - len(matched_b2)

            print(
                f"Matched pairs: {len(matches)} | "
                f"Unmatched B1: {unmatched_b1} | "
                f"Unmatched B2: {unmatched_b2} | "
                f"Mean IoU: {mean_iou:.4f}"
            )

            for match in matches:
                match_rows.append({
                    "date": date,
                    "threshold_multiplier": multiplier,
                    **match,
                })

            # One summary row per date and multiplier.
            summary_rows.append({
                "date": date,
                "threshold_multiplier": multiplier,

                "band1_threshold": detected[1]["threshold"],
                "band1_robust_sigma": detected[1]["robust_sigma"],
                "band1_candidate_samples": detected[1]["candidate_count"],
                "band1_flagged_percent": detected[1]["candidate_percent"],
                "band1_events": len(band1_events),

                "band2_threshold": detected[2]["threshold"],
                "band2_robust_sigma": detected[2]["robust_sigma"],
                "band2_candidate_samples": detected[2]["candidate_count"],
                "band2_flagged_percent": detected[2]["candidate_percent"],
                "band2_events": len(band2_events),

                "matched_pairs": len(matches),
                "band1_unmatched": unmatched_b1,
                "band2_unmatched": unmatched_b2,
                "mean_iou": mean_iou,
            })

    # ========================================================
    # SAVE RESULTS
    # ========================================================

    summary_df = pd.DataFrame(summary_rows)
    events_df = pd.DataFrame(event_rows)
    matches_df = pd.DataFrame(match_rows)

    summary_path = (
        OUTPUT_DIR / "threshold_sensitivity_summary.csv"
    )

    events_path = (
        OUTPUT_DIR / "threshold_sensitivity_events.csv"
    )

    matches_path = (
        OUTPUT_DIR / "threshold_sensitivity_matches.csv"
    )

    summary_df.to_csv(summary_path, index=False)
    events_df.to_csv(events_path, index=False)
    matches_df.to_csv(matches_path, index=False)

    print("\n" + "=" * 65)
    print("THRESHOLD SENSITIVITY EXPERIMENT COMPLETE")
    print("=" * 65)

    print("\nSummary:")
    print(
        summary_df[
            [
                "date",
                "threshold_multiplier",
                "band1_flagged_percent",
                "band1_events",
                "band2_flagged_percent",
                "band2_events",
                "matched_pairs",
                "mean_iou",
            ]
        ].to_string(index=False)
    )

    print("\nSaved:")
    print(summary_path)
    print(events_path)
    print(matches_path)


if __name__ == "__main__":
    main()
