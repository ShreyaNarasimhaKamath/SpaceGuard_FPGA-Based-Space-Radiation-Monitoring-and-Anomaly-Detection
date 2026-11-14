
from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset


# --------------------------------------------------
# CONFIGURATION
# --------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = ROOT / "data" / "raw" / "goes16_mpshi"
OUTPUT_DIR = ROOT / "results" / "evaluation"

DATES = ["20240101", "20240102", "20240103"]

TELESCOPE_INDEX = 0
BANDS = {
    1: 0,
    2: 1,
}

WINDOW = 1800
CALIBRATION_SECONDS = 10 * 3600

THRESHOLD_MULTIPLIER = 6.0
MERGE_GAP_SECONDS = 120
MIN_ACTIVE_SAMPLES = 60

# EMA uses the same nominal 1800-sample time scale.
EMA_ALPHA = 2.0 / (WINDOW + 1)


# --------------------------------------------------
# DATA LOADING
# --------------------------------------------------

def load_channel(date, band_index):

    filename = DATA_DIR / f"goes16_mpshi_{date}.nc"

    if not filename.exists():
        raise FileNotFoundError(filename)

    with Dataset(filename, "r") as ds:

        flux_raw = ds.variables["DiffElectronFluxes"][
            :, TELESCOPE_INDEX, band_index
        ]

        dqf_raw = ds.variables["DiffElectronFluxDQFs"][
            :, TELESCOPE_INDEX, band_index
        ]

        # Convert masked arrays to floating point before filling.
        flux = np.ma.asarray(flux_raw).astype(float).filled(np.nan)
        dqf = np.ma.asarray(dqf_raw).astype(float).filled(np.nan)

    valid = (
        np.isfinite(flux)
        & np.isfinite(dqf)
        & (dqf == 0)
        & (flux >= 0)
    )

    signal = pd.Series(flux)
    signal.loc[~valid] = np.nan

    return signal, valid


# --------------------------------------------------
# BASELINE ESTIMATORS
# --------------------------------------------------

def calculate_baselines(signal):

    rolling_median = signal.rolling(
        window=WINDOW,
        min_periods=WINDOW
    ).median()

    moving_average = signal.rolling(
        window=WINDOW,
        min_periods=WINDOW
    ).mean()

    ema_always = signal.ewm(
        alpha=EMA_ALPHA,
        adjust=False,
        min_periods=WINDOW,
        ignore_na=True
    ).mean()

    return {
        "rolling_median": rolling_median,
        "moving_average": moving_average,
        "ema_always": ema_always,
    }


# --------------------------------------------------
# ROBUST THRESHOLD
# --------------------------------------------------

def calculate_threshold(residual, calibration_mask):

    calibration = residual[calibration_mask]
    calibration = calibration[np.isfinite(calibration)]

    if len(calibration) < 100:
        raise ValueError(
            f"Insufficient common calibration samples: {len(calibration)}"
        )

    median = np.median(calibration)
    mad = np.median(np.abs(calibration - median))

    robust_sigma = 1.4826 * mad
    threshold = THRESHOLD_MULTIPLIER * robust_sigma

    return median, mad, robust_sigma, threshold, len(calibration)


# --------------------------------------------------
# EVENT GROUPING
# --------------------------------------------------

def count_events(candidate_mask):

    indices = np.flatnonzero(candidate_mask)

    if len(indices) == 0:
        return 0, []

    events = []

    start = indices[0]
    previous = indices[0]
    active_count = 1

    for idx in indices[1:]:

        gap = idx - previous - 1

        if gap <= MERGE_GAP_SECONDS:
            active_count += 1
        else:
            duration = previous - start + 1

            if active_count >= MIN_ACTIVE_SAMPLES:
                events.append({
                    "start_sample": int(start),
                    "end_sample": int(previous),
                    "duration_seconds": int(duration),
                    "active_samples": int(active_count),
                })

            start = idx
            active_count = 1

        previous = idx

    duration = previous - start + 1

    if active_count >= MIN_ACTIVE_SAMPLES:
        events.append({
            "start_sample": int(start),
            "end_sample": int(previous),
            "duration_seconds": int(duration),
            "active_samples": int(active_count),
        })

    return len(events), events


# --------------------------------------------------
# MAIN COMPARISON
# --------------------------------------------------

def main():

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    event_rows = []

    methods = [
        "rolling_median",
        "moving_average",
        "ema_always",
    ]

    for date in DATES:

        print(f"\nProcessing {date}")

        for band, band_index in BANDS.items():

            print(f"  Band {band}")

            signal, input_valid = load_channel(date, band_index)

            baselines = calculate_baselines(signal)

            residuals = {
                method: signal - baseline
                for method, baseline in baselines.items()
            }

            # --------------------------------------------------
            # COMMON EVALUATION MASK
            # --------------------------------------------------

            common_mask = input_valid.copy()

            for method in methods:
                common_mask &= np.isfinite(
                    baselines[method].to_numpy()
                )

                common_mask &= np.isfinite(
                    residuals[method].to_numpy()
                )

            common_indices = np.flatnonzero(common_mask)

            if len(common_indices) == 0:
                raise ValueError(
                    f"No common evaluation samples: {date}, Band {band}"
                )

            # First 10 hours, restricted to common timestamps.
            calibration_mask = (
                common_mask
                & (np.arange(len(signal)) < CALIBRATION_SECONDS)
            )

            # Use the same calibration timestamps for every estimator.
            reference_calibration = residuals[
                "rolling_median"
            ].to_numpy()[calibration_mask]

            print(
                f"    Common valid samples: {len(common_indices)}"
            )
            print(
                f"    Common calibration samples: "
                f"{calibration_mask.sum()}"
            )

            # --------------------------------------------------
            # CALCULATE DETECTIONS
            # --------------------------------------------------

            candidate_masks = {}
            thresholds = {}

            for method in methods:

                residual = residuals[method].to_numpy()

                median, mad, sigma, threshold, cal_count = (
                    calculate_threshold(
                        residual,
                        calibration_mask
                    )
                )

                # Positive excursions only, consistent with reference.
                candidate = (
                    common_mask
                    & (residual > threshold)
                )

                candidate_masks[method] = candidate
                thresholds[method] = threshold

            reference_mask = candidate_masks["rolling_median"]

            # --------------------------------------------------
            # METRICS
            # --------------------------------------------------

            for method in methods:

                candidate = candidate_masks[method]

                intersection = np.sum(
                    candidate & reference_mask
                )

                union = np.sum(
                    candidate | reference_mask
                )

                iou = (
                    intersection / union
                    if union > 0
                    else 1.0
                )

                candidate_count = int(candidate.sum())

                flagged_percentage = (
                    100.0 * candidate_count / len(common_indices)
                )

                event_count, events = count_events(candidate)

                residual = residuals[method].to_numpy()

                median, mad, sigma, threshold, cal_count = (
                    calculate_threshold(
                        residual,
                        calibration_mask
                    )
                )

                summary_rows.append({
                    "date": date,
                    "band": band,
                    "method": method,
                    "input_valid_samples": int(input_valid.sum()),
                    "common_evaluation_samples": len(common_indices),
                    "common_calibration_samples": cal_count,
                    "baseline_available_samples": int(
                        np.isfinite(
                            baselines[method].to_numpy()
                        ).sum()
                    ),
                    "calibration_median": median,
                    "calibration_mad": mad,
                    "robust_sigma": sigma,
                    "threshold": threshold,
                    "candidate_samples": candidate_count,
                    "flagged_percentage": flagged_percentage,
                    "event_count": event_count,
                    "sample_mask_iou_vs_median": iou,
                })

                for event_number, event in enumerate(events, start=1):

                    event_rows.append({
                        "date": date,
                        "band": band,
                        "method": method,
                        "event_number": event_number,
                        **event,
                    })

                print(
                    f"    {method:16s} "
                    f"flagged={flagged_percentage:6.2f}% "
                    f"events={event_count:3d} "
                    f"IoU={iou:.3f}"
                )

    # --------------------------------------------------
    # SAVE OUTPUTS
    # --------------------------------------------------

    summary = pd.DataFrame(summary_rows)
    events = pd.DataFrame(event_rows)

    summary_path = OUTPUT_DIR / "fair_baseline_comparison.csv"
    events_path = OUTPUT_DIR / "fair_event_comparison.csv"

    summary.to_csv(summary_path, index=False)
    events.to_csv(events_path, index=False)

    print("\nFAIR BASELINE COMPARISON COMPLETE")
    print(f"\nSaved: {summary_path}")
    print(f"Saved: {events_path}")

    print("\nSUMMARY")
    print(
        summary[
            [
                "date",
                "band",
                "method",
                "common_evaluation_samples",
                "common_calibration_samples",
                "flagged_percentage",
                "event_count",
                "sample_mask_iou_vs_median",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
