
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset


# ============================================================
# CONFIGURATION
# ============================================================

ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = ROOT / "data" / "raw" / "goes16_mpshi"
PLOT_DIR = ROOT / "results" / "plots"
EVAL_DIR = ROOT / "results" / "evaluation"

PLOT_DIR.mkdir(parents=True, exist_ok=True)
EVAL_DIR.mkdir(parents=True, exist_ok=True)

DATES = ["20240101", "20240102", "20240103"]

TELESCOPE_INDEX = 0
BANDS = [0, 1]              # Band 1 and Band 2

WINDOW = 1800               # 30 minutes at 1 Hz
CALIBRATION_SECONDS = 36000  # First 10 hours

THRESHOLD_MULTIPLIER = 6.0
EMA_ALPHA = 2.0 / (WINDOW + 1)

MERGE_GAP_SECONDS = 120
MIN_ACTIVE_SAMPLES = 60


# ============================================================
# DATA LOADING
# ============================================================

def load_channel(date, band_index):

    filename = DATA_DIR / f"goes16_mpshi_{date}.nc"

    if not filename.exists():
        raise FileNotFoundError(f"Missing dataset: {filename}")

    with Dataset(filename, "r") as nc:

        flux_raw = nc.variables["DiffElectronFluxes"][:]
        dqf_raw = nc.variables["DiffElectronFluxDQFs"][:]

        flux = np.asarray(
            np.ma.filled(flux_raw, np.nan),
            dtype=np.float64
        )[:, TELESCOPE_INDEX, band_index]

        dqf = np.asarray(
    np.ma.filled(dqf_raw.astype(np.int16), -1),
    dtype=np.int16
)[:, TELESCOPE_INDEX, band_index]

    valid = (
        (dqf == 0)
        & np.isfinite(flux)
        & (flux >= 0)
    )

    return flux, valid


# ============================================================
# ROBUST THRESHOLD
# ============================================================

def calculate_threshold(residual, valid, date):

    calibration_mask = (
        (np.arange(len(residual)) < CALIBRATION_SECONDS)
        & valid
        & np.isfinite(residual)
    )

    calibration = residual[calibration_mask]

    if len(calibration) == 0:
        raise ValueError(f"No calibration residuals for {date}")

    center = np.median(calibration)
    mad = np.median(np.abs(calibration - center))

    robust_sigma = 1.4826 * mad
    threshold = THRESHOLD_MULTIPLIER * robust_sigma

    return center, mad, robust_sigma, threshold, len(calibration)


# ============================================================
# ROLLING BASELINES
# ============================================================

def rolling_baseline(signal, valid, method):

    series = pd.Series(
        np.where(valid, signal, np.nan)
    )

    if method == "rolling_median":
        baseline = series.rolling(
            window=WINDOW,
            min_periods=WINDOW
        ).median()

    elif method == "moving_average":
        baseline = series.rolling(
            window=WINDOW,
            min_periods=WINDOW
        ).mean()

    else:
        raise ValueError(f"Unknown rolling method: {method}")

    return baseline.to_numpy(dtype=np.float64)


# ============================================================
# EMA BASELINE
# ============================================================

def ema_baseline(signal, valid, freeze_threshold=None):

    """
    Causal EMA.

    The baseline at sample n is the state formed from previous
    samples. The current sample updates the state after its
    residual has been calculated.

    If freeze_threshold is supplied, the state is frozen when
    a positive residual exceeds that threshold.

    The first 1800 valid samples are used for warm-up.
    """

    baseline = np.full(len(signal), np.nan)

    state = None
    valid_count = 0

    for i in range(len(signal)):

        if not valid[i]:
            continue

        x = signal[i]

        if state is None:
            state = x

        valid_count += 1

        # Baseline becomes available after 1800 valid samples.
        if valid_count >= WINDOW:
            baseline[i] = state

        # During the calibration period, always update.
        # After calibration, optionally freeze on candidates.
        freeze = False

        if (
            freeze_threshold is not None
            and i >= CALIBRATION_SECONDS
            and valid_count >= WINDOW
        ):
            residual_before_update = x - state
            freeze = residual_before_update > freeze_threshold

        if not freeze:
            state = state + EMA_ALPHA * (x - state)

    return baseline


# ============================================================
# EVENT GROUPING
# ============================================================

def count_events(candidate):

    indices = np.flatnonzero(candidate)

    if len(indices) == 0:
        return 0

    groups = []

    start = indices[0]
    previous = indices[0]
    active_count = 1

    for index in indices[1:]:

        gap = index - previous - 1

        if gap <= MERGE_GAP_SECONDS:
            active_count += 1

        else:
            groups.append(
                (start, previous, active_count)
            )

            start = index
            active_count = 1

        previous = index

    groups.append((start, previous, active_count))

    # Retain only events with sufficient active samples.
    accepted = [
        event for event in groups
        if event[2] >= MIN_ACTIVE_SAMPLES
    ]

    return len(accepted)


# ============================================================
# DETECTOR METRICS
# ============================================================

def evaluate_detector(
    date,
    band_number,
    method,
    signal,
    valid,
    baseline
):

    residual = signal - baseline

    usable = (
        valid
        & np.isfinite(baseline)
        & np.isfinite(residual)
    )

    center, mad, robust_sigma, threshold, cal_count = (
        calculate_threshold(residual, usable, date)
    )

    candidate = usable & (residual > threshold)

    usable_count = int(np.sum(usable))
    candidate_count = int(np.sum(candidate))

    flagged_percentage = (
        100.0 * candidate_count / usable_count
        if usable_count else 0.0
    )

    event_count = count_events(candidate)

    return {
        "date": date,
        "band": band_number,
        "method": method,
        "valid_samples": int(np.sum(valid)),
        "baseline_available_samples": usable_count,
        "calibration_samples": cal_count,
        "calibration_median": center,
        "calibration_mad": mad,
        "robust_sigma": robust_sigma,
        "threshold": threshold,
        "candidate_samples": candidate_count,
        "flagged_percentage": flagged_percentage,
        "event_count": event_count,
    }, residual, candidate


# ============================================================
# MAIN EXPERIMENT
# ============================================================

def main():

    all_metrics = []
    all_results = {}

    for date in DATES:

        print(f"\nProcessing {date}")

        for band_index in BANDS:

            band_number = band_index + 1

            print(f"  Band {band_number}")

            signal, valid = load_channel(date, band_index)

            # ------------------------------------------------
            # 1. Reference rolling median
            # ------------------------------------------------

            median_base = rolling_baseline(
                signal, valid, "rolling_median"
            )

            median_metrics, median_residual, median_candidate = (
                evaluate_detector(
                    date,
                    band_number,
                    "rolling_median",
                    signal,
                    valid,
                    median_base
                )
            )

            # ------------------------------------------------
            # 2. Causal moving average
            # ------------------------------------------------

            mean_base = rolling_baseline(
                signal, valid, "moving_average"
            )

            mean_metrics, mean_residual, mean_candidate = (
                evaluate_detector(
                    date,
                    band_number,
                    "moving_average",
                    signal,
                    valid,
                    mean_base
                )
            )

            # ------------------------------------------------
            # 3. EMA, always updating
            # ------------------------------------------------

            ema_base = ema_baseline(signal, valid)

            ema_metrics, ema_residual, ema_candidate = (
                evaluate_detector(
                    date,
                    band_number,
                    "ema_always",
                    signal,
                    valid,
                    ema_base
                )
            )

            # ------------------------------------------------
            # 4. EMA, freeze when candidate is detected
            # ------------------------------------------------

            # First obtain the EMA's own calibration threshold.
            ema_threshold = ema_metrics["threshold"]

            frozen_ema_base = ema_baseline(
                signal,
                valid,
                freeze_threshold=ema_threshold
            )

            frozen_metrics, frozen_residual, frozen_candidate = (
                evaluate_detector(
                    date,
                    band_number,
                    "ema_freeze",
                    signal,
                    valid,
                    frozen_ema_base
                )
            )

            # ------------------------------------------------
            # Compare candidate masks with median reference
            # ------------------------------------------------

            configurations = [
                (
                    median_metrics,
                    median_base,
                    median_residual,
                    median_candidate
                ),
                (
                    mean_metrics,
                    mean_base,
                    mean_residual,
                    mean_candidate
                ),
                (
                    ema_metrics,
                    ema_base,
                    ema_residual,
                    ema_candidate
                ),
                (
                    frozen_metrics,
                    frozen_ema_base,
                    frozen_residual,
                    frozen_candidate
                )
            ]

            reference = median_candidate

            for metrics, baseline, residual, candidate in configurations:

                intersection = np.sum(candidate & reference)
                union = np.sum(candidate | reference)

                jaccard = (
                    intersection / union
                    if union > 0 else 1.0
                )

                metrics["sample_mask_iou_vs_median"] = jaccard

                all_metrics.append(metrics)

                all_results[(date, band_number, metrics["method"])] = (
                    signal,
                    baseline,
                    residual,
                    candidate
                )

                print(
                    f"    {metrics['method']:16s} "
                    f"flagged={metrics['flagged_percentage']:6.2f}% "
                    f"events={metrics['event_count']:3d} "
                    f"IoU={jaccard:.3f}"
                )

    # ========================================================
    # SAVE METRICS
    # ========================================================

    metrics_df = pd.DataFrame(all_metrics)

    metrics_path = (
        EVAL_DIR / "baseline_estimator_comparison.csv"
    )

    metrics_df.to_csv(metrics_path, index=False)

    print("\nSaved metrics:")
    print(metrics_path)

    # ========================================================
    # PLOT 1: BASELINE COMPARISON
    # ========================================================

    methods = [
        "rolling_median",
        "moving_average",
        "ema_always",
        "ema_freeze"
    ]

    fig, axes = plt.subplots(
        len(DATES),
        len(BANDS),
        figsize=(17, 13),
        sharex=False
    )

    for row, date in enumerate(DATES):

        for col, band_index in enumerate(BANDS):

            band_number = band_index + 1
            ax = axes[row, col]

            signal = all_results[
                (date, band_number, "rolling_median")
            ][0]

            hours = np.arange(len(signal)) / 3600.0

            # Plot raw signal.
            ax.plot(
                hours,
                signal,
                color="black",
                linewidth=0.7,
                alpha=0.5,
                label="Original flux"
            )

            for method in methods:

                baseline = all_results[
                    (date, band_number, method)
                ][1]

                ax.plot(
                    hours,
                    baseline,
                    linewidth=1.0,
                    label=method
                )

            ax.axvline(
                CALIBRATION_SECONDS / 3600,
                linestyle=":",
                color="red",
                label="Calibration boundary"
            )

            ax.set_title(
                f"{date} | Band {band_number}"
            )

            ax.set_xlabel("Hour of day")
            ax.set_ylabel("Flux")
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=7)

    fig.suptitle(
        "GOES-16 Baseline Estimator Comparison",
        fontsize=15
    )

    fig.tight_layout()

    plot_path = (
        PLOT_DIR / "baseline_estimator_comparison.png"
    )

    fig.savefig(plot_path, dpi=180)
    plt.close(fig)

    # ========================================================
    # PLOT 2: FLAGGED SAMPLE PERCENTAGE
    # ========================================================

    fig, axes = plt.subplots(
        len(DATES),
        len(BANDS),
        figsize=(14, 11)
    )

    for row, date in enumerate(DATES):

        for col, band_index in enumerate(BANDS):

            band_number = band_index + 1

            subset = metrics_df[
                (metrics_df["date"] == date)
                & (metrics_df["band"] == band_number)
            ]

            ax = axes[row, col]

            ax.bar(
                subset["method"],
                subset["flagged_percentage"]
            )

            ax.set_title(
                f"{date} | Band {band_number}"
            )

            ax.set_ylabel("Flagged samples (%)")
            ax.tick_params(axis="x", rotation=25)
            ax.grid(axis="y", alpha=0.3)

    fig.suptitle(
        "Candidate Sample Percentage by Baseline Estimator"
    )

    fig.tight_layout()

    percentage_path = (
        PLOT_DIR / "baseline_estimator_flagged_percentage.png"
    )

    fig.savefig(percentage_path, dpi=180)
    plt.close(fig)

    print("\nGenerated plots:")
    print(plot_path)
    print(percentage_path)

    print("\nBASELINE COMPARISON COMPLETE")


if __name__ == "__main__":
    main()
