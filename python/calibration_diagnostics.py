
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset


# ============================================================
# CONFIGURATION — MATCH THE REFERENCE DETECTOR
# ============================================================

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

BASELINE_WINDOW = 1800
CALIBRATION_HOURS = 10
CALIBRATION_SAMPLES = CALIBRATION_HOURS * 3600

REFERENCE_MULTIPLIER = 6
TEST_MULTIPLIERS = [4, 5, 6, 7, 8]

MERGE_GAP = 120
MIN_ACTIVE_SAMPLES = 60

PLOT_DIR.mkdir(parents=True, exist_ok=True)
EVAL_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# DATA LOADING
# ============================================================

def load_channel(path, band_index):

    with Dataset(path, "r") as nc:

        flux_data = nc.variables["DiffElectronFluxes"][
            :, 0, band_index
        ]

        dqf_data = nc.variables["DiffElectronFluxDQFs"][
            :, 0, band_index
        ]

        flux = (
            np.ma.asarray(flux_data)
            .astype(float)
            .filled(np.nan)
        )

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

    return flux, dqf, valid


# ============================================================
# RESIDUAL AND THRESHOLD CALCULATION
# ============================================================

def calculate_residual(flux, valid):

    signal = pd.Series(
        np.where(valid, flux, np.nan)
    )

    baseline = signal.rolling(
        window=BASELINE_WINDOW,
        min_periods=BASELINE_WINDOW
    ).median()

    residual = signal - baseline

    calibration = residual.iloc[
        :CALIBRATION_SAMPLES
    ].dropna()

    if calibration.empty:
        raise ValueError("No valid calibration residual samples.")

    calibration_median = calibration.median()

    calibration_mad = np.median(
        np.abs(
            calibration.to_numpy() - calibration_median
        )
    )

    robust_sigma = 1.4826 * calibration_mad

    thresholds = {
        multiplier: multiplier * robust_sigma
        for multiplier in TEST_MULTIPLIERS
    }

    return (
        residual.to_numpy(),
        calibration.to_numpy(),
        robust_sigma,
        thresholds
    )


# ============================================================
# EVENT GROUPING
# ============================================================

def group_events(candidate):

    indices = np.flatnonzero(candidate)

    if len(indices) == 0:
        return []

    groups = []

    start = indices[0]
    previous = indices[0]
    active_count = 1

    for index in indices[1:]:

        gap = index - previous - 1

        if gap <= MERGE_GAP:
            active_count += 1

        else:
            groups.append(
                (start, previous, active_count)
            )

            start = index
            active_count = 1

        previous = index

    groups.append(
        (start, previous, active_count)
    )

    return [
        {
            "start_index": int(start),
            "end_index": int(end),
            "active_samples": int(active),
        }
        for start, end, active in groups
        if active >= MIN_ACTIVE_SAMPLES
    ]


# ============================================================
# INVALID SAMPLE INTERVALS
# ============================================================

def extract_invalid_intervals(valid, dqf):

    invalid = ~valid
    indices = np.flatnonzero(invalid)

    if len(indices) == 0:
        return []

    intervals = []

    start = indices[0]
    previous = indices[0]

    for index in indices[1:]:

        if index != previous + 1:

            interval = {
                "start_index": int(start),
                "end_index": int(previous),
                "start_hour": start / 3600,
                "end_hour": previous / 3600,
                "sample_count": int(previous - start + 1),
            }

            interval_dqf = dqf[start:previous + 1]

            values, counts = np.unique(
                interval_dqf,
                return_counts=True
            )

            interval["dqf_values"] = ",".join(
                str(v) for v in values
            )

            intervals.append(interval)

            start = index

        previous = index

    interval = {
        "start_index": int(start),
        "end_index": int(previous),
        "start_hour": start / 3600,
        "end_hour": previous / 3600,
        "sample_count": int(previous - start + 1),
    }

    interval_dqf = dqf[start:previous + 1]

    values, counts = np.unique(
        interval_dqf,
        return_counts=True
    )

    interval["dqf_values"] = ",".join(
        str(v) for v in values
    )

    intervals.append(interval)

    return intervals


# ============================================================
# MAIN ANALYSIS
# ============================================================

def main():

    residual_records = []
    invalid_records = []
    interval_records = []
    event_margin_records = []

    plot_data = {}

    for filename in FILES:

        date = filename.split("_")[-1].split(".")[0]
        path = DATA_DIR / filename

        if not path.exists():
            raise FileNotFoundError(path)

        print("\n" + "=" * 65)
        print(f"PROCESSING {date}")
        print("=" * 65)

        for band, band_index in BANDS.items():

            flux, dqf, valid = load_channel(
                path,
                band_index
            )

            residual, calibration, sigma, thresholds = (
                calculate_residual(flux, valid)
            )

            calibration_residual = residual[
                :CALIBRATION_SAMPLES
            ]

            post_calibration_residual = residual[
                CALIBRATION_SAMPLES:
            ]

            calibration_valid = (
                np.isfinite(calibration_residual)
            )

            post_valid = (
                np.isfinite(post_calibration_residual)
            )

            calibration_values = calibration_residual[
                calibration_valid
            ]

            post_values = post_calibration_residual[
                post_valid
            ]

            threshold = thresholds[
                REFERENCE_MULTIPLIER
            ]

            # ------------------------------------------------
            # Residual statistics
            # ------------------------------------------------

            residual_records.append({

                "date": date,
                "band": band,

                "valid_flux_samples": int(valid.sum()),
                "invalid_flux_samples": int((~valid).sum()),

                "calibration_residual_samples": len(
                    calibration_values
                ),

                "post_calibration_residual_samples": len(
                    post_values
                ),

                "calibration_median": float(
                    np.median(calibration_values)
                ),

                "calibration_std": float(
                    np.std(calibration_values)
                ),

                "calibration_mad": float(
                    np.median(
                        np.abs(
                            calibration_values
                            - np.median(calibration_values)
                        )
                    )
                ),

                "post_calibration_median": float(
                    np.median(post_values)
                ),

                "post_calibration_std": float(
                    np.std(post_values)
                ),

                "post_calibration_95th": float(
                    np.percentile(post_values, 95)
                ),

                "post_calibration_99th": float(
                    np.percentile(post_values, 99)
                ),

                "robust_sigma": float(sigma),
                "threshold_6x": float(threshold),

                "calibration_samples_above_6x": int(
                    np.sum(calibration_values > threshold)
                ),

                "post_calibration_samples_above_6x": int(
                    np.sum(post_values > threshold)
                ),
            })

            plot_data[(date, band)] = {
                "calibration": calibration_values,
                "post": post_values,
                "threshold": threshold,
            }

            # ------------------------------------------------
            # Invalid sample statistics
            # ------------------------------------------------

            intervals = extract_invalid_intervals(
                valid,
                dqf
            )

            invalid_records.append({

                "date": date,
                "band": band,
                "total_samples": len(valid),
                "valid_samples": int(valid.sum()),
                "invalid_samples": int((~valid).sum()),

                "invalid_percent": float(
                    100 * (~valid).sum() / len(valid)
                ),

                "invalid_intervals": len(intervals),
                "first_invalid_hour": (
                    intervals[0]["start_hour"]
                    if intervals else np.nan
                ),
            })

            for interval_number, interval in enumerate(
                intervals, start=1
            ):

                interval_records.append({
                    "date": date,
                    "band": band,
                    "interval_id": interval_number,
                    **interval,
                })

            # ------------------------------------------------
            # Threshold-margin analysis
            # ------------------------------------------------

            candidates_6x = (
                np.isfinite(residual)
                & valid
                & (residual > threshold)
            )

            reference_events = group_events(
                candidates_6x
            )

            for event_id, event in enumerate(
                reference_events, start=1
            ):

                start = event["start_index"]
                end = event["end_index"]

                # Only valid residual samples inside event span.
                event_residual = residual[start:end + 1]
                event_valid = valid[start:end + 1]

                usable = (
                    np.isfinite(event_residual)
                    & event_valid
                )

                values = event_residual[usable]

                if len(values) == 0:
                    continue

                row = {

                    "date": date,
                    "band": band,
                    "reference_event_id": event_id,

                    "start_hour": start / 3600,
                    "end_hour": end / 3600,

                    "active_samples": event["active_samples"],

                    "peak_residual": float(np.max(values)),
                    "median_residual": float(np.median(values)),

                    "threshold_6x": float(threshold),

                    "peak_to_threshold_ratio": float(
                        np.max(values) / threshold
                    ),
                }

                for multiplier in TEST_MULTIPLIERS:

                    current_threshold = thresholds[multiplier]

                    row[
                        f"fraction_above_{multiplier}x"
                    ] = float(
                        np.mean(values > current_threshold)
                    )

                event_margin_records.append(row)

            print(
                f"Band {band}: "
                f"valid={valid.sum()}, "
                f"invalid={(~valid).sum()}, "
                f"invalid intervals={len(intervals)}, "
                f"calibration residuals={len(calibration_values)}, "
                f"post-calibration residuals={len(post_values)}, "
                f"6x events={len(reference_events)}"
            )

    # ========================================================
    # PLOT 1: CALIBRATION VS POST-CALIBRATION DISTRIBUTIONS
    # ========================================================

    fig, axes = plt.subplots(
        3, 2,
        figsize=(14, 13),
        constrained_layout=True
    )

    for row, filename in enumerate(FILES):

        date = filename.split("_")[-1].split(".")[0]

        for col, band in enumerate([1, 2]):

            ax = axes[row, col]

            data = plot_data[(date, band)]

            ax.hist(
                data["calibration"],
                bins=120,
                density=True,
                alpha=0.55,
                label="First 10 hours"
            )

            ax.hist(
                data["post"],
                bins=120,
                density=True,
                alpha=0.45,
                label="Remaining 14 hours"
            )

            ax.axvline(
                data["threshold"],
                linestyle="--",
                label="Positive 6x threshold"
            )

            ax.set_title(f"{date} | Band {band}")
            ax.set_xlabel("Residual flux")
            ax.set_ylabel("Density")
            ax.legend(fontsize=8)

    fig.suptitle(
        "Calibration vs Post-Calibration Residual Distributions"
    )

    fig.savefig(
        PLOT_DIR / "calibration_vs_postcalibration.png",
        dpi=200
    )

    plt.close(fig)

    # ========================================================
    # PLOT 2: INVALID SAMPLE TIMELINE
    # ========================================================

    fig, axes = plt.subplots(
        3, 2,
        figsize=(14, 8),
        constrained_layout=True
    )

    for row, filename in enumerate(FILES):

        date = filename.split("_")[-1].split(".")[0]

        for col, band in enumerate([1, 2]):

            ax = axes[row, col]

            matching = [
                record for record in interval_records
                if record["date"] == date
                and record["band"] == band
            ]

            for interval in matching:

                ax.axvspan(
                    interval["start_hour"],
                    (interval["end_index"] + 1) / 3600,
                    alpha=0.65
                )

            ax.set_xlim(0, 24)
            ax.set_xticks(range(0, 25, 2))
            ax.set_title(
                f"{date} | Band {band} invalid intervals"
            )
            ax.set_xlabel("Hour of day")
            ax.set_ylabel("Invalid sample periods")

    fig.suptitle("Invalid Sample Timeline")

    fig.savefig(
        PLOT_DIR / "invalid_sample_timeline.png",
        dpi=200
    )

    plt.close(fig)

    # ========================================================
    # SAVE CSV FILES
    # ========================================================

    residual_df = pd.DataFrame(residual_records)
    invalid_df = pd.DataFrame(invalid_records)
    interval_df = pd.DataFrame(interval_records)
    margin_df = pd.DataFrame(event_margin_records)

    residual_df.to_csv(
        EVAL_DIR / "calibration_statistics.csv",
        index=False
    )

    invalid_df.to_csv(
        EVAL_DIR / "invalid_sample_summary.csv",
        index=False
    )

    interval_df.to_csv(
        EVAL_DIR / "invalid_sample_intervals.csv",
        index=False
    )

    margin_df.to_csv(
        EVAL_DIR / "event_threshold_margins.csv",
        index=False
    )

    # ========================================================
    # PRINT RESULTS
    # ========================================================

    print("\n" + "=" * 65)
    print("CALIBRATION DIAGNOSTICS COMPLETE")
    print("=" * 65)

    print("\nInvalid sample summary:")
    print(invalid_df.to_string(index=False))

    print("\nCalibration statistics:")
    print(
        residual_df[
            [
                "date",
                "band",
                "calibration_residual_samples",
                "post_calibration_residual_samples",
                "calibration_median",
                "post_calibration_median",
                "robust_sigma",
                "threshold_6x",
            ]
        ].to_string(index=False)
    )

    print("\nGenerated plots:")
    print(PLOT_DIR / "calibration_vs_postcalibration.png")
    print(PLOT_DIR / "invalid_sample_timeline.png")

    print("\nGenerated evaluation files:")
    print(EVAL_DIR / "calibration_statistics.csv")
    print(EVAL_DIR / "invalid_sample_summary.csv")
    print(EVAL_DIR / "invalid_sample_intervals.csv")
    print(EVAL_DIR / "event_threshold_margins.csv")


if __name__ == "__main__":
    main()
