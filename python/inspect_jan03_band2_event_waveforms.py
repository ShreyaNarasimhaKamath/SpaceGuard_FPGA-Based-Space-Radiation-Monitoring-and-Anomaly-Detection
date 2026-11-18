
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset


# --------------------------------------------------
# CONFIGURATION
# --------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]

DATA_FILE = (
    ROOT / "data" / "raw" / "goes16_mpshi"
    / "goes16_mpshi_20240103.nc"
)

PLOT_DIR = ROOT / "results" / "plots"
RESULT_DIR = ROOT / "results" / "evaluation"

TELESCOPE_INDEX = 0
BAND_INDEX = 1  # Electron Band 2

WINDOWS = [300, 600, 900, 1800]
REFERENCE_WINDOW = 1800

CALIBRATION_SAMPLES = 10 * 3600
THRESHOLD_MULTIPLIER = 6

# Reference rolling-median events from Jan 3, Band 2.
REFERENCE_EVENTS = [
    (20525, 21160),
    (37459, 38088),
    (38513, 38692),
    (57039, 57709),
    (73603, 74213),
]

# Context shown before and after each event.
CONTEXT_SECONDS = 900


# --------------------------------------------------
# LOAD DATA
# --------------------------------------------------

def load_data():
    if not DATA_FILE.exists():
        raise FileNotFoundError(DATA_FILE)

    with Dataset(DATA_FILE, "r") as nc:
        flux_raw = nc.variables["DiffElectronFluxes"][:]
        dqf_raw = nc.variables["DiffElectronFluxDQFs"][:]

    flux = np.ma.asarray(flux_raw).astype(float).filled(np.nan)
    dqf = np.ma.asarray(dqf_raw).astype(float).filled(np.nan)

    signal = flux[:, TELESCOPE_INDEX, BAND_INDEX]
    dqf_channel = dqf[:, TELESCOPE_INDEX, BAND_INDEX]

    valid = (
        np.isfinite(signal)
        & np.isfinite(dqf_channel)
        & (dqf_channel == 0)
        & (signal >= 0)
    )

    # Jan 3 Band 2 should be fully valid, but verify.
    print(f"Total samples: {len(signal)}")
    print(f"Valid samples: {valid.sum()}")
    print(f"Invalid samples: {(~valid).sum()}")

    return signal, valid


# --------------------------------------------------
# BASELINES
# --------------------------------------------------

def calculate_baselines(signal):
    series = pd.Series(signal)

    baselines = {
        "Rolling median (30 min)": series.rolling(
            window=REFERENCE_WINDOW,
            min_periods=REFERENCE_WINDOW
        ).median().to_numpy()
    }

    for window in WINDOWS:
        baselines[f"Moving average ({window // 60} min)"] = (
            series.rolling(
                window=window,
                min_periods=window
            ).mean().to_numpy()
        )

    return baselines


# --------------------------------------------------
# CALIBRATION THRESHOLDS
# --------------------------------------------------

def calculate_thresholds(signal, valid, baselines):
    thresholds = {}
    calibration_indices = np.arange(len(signal)) < CALIBRATION_SAMPLES

    for name, baseline in baselines.items():
        residual = signal - baseline

        calibration_mask = (
            valid
            & calibration_indices
            & np.isfinite(residual)
        )

        calibration = residual[calibration_mask]

        median = np.median(calibration)
        mad = np.median(np.abs(calibration - median))
        robust_sigma = 1.4826 * mad

        threshold = THRESHOLD_MULTIPLIER * robust_sigma

        thresholds[name] = {
            "median": median,
            "mad": mad,
            "robust_sigma": robust_sigma,
            "threshold": threshold,
        }

    return thresholds


# --------------------------------------------------
# PLOT EVENT WAVEFORMS
# --------------------------------------------------

def plot_event(
    event_number,
    start,
    end,
    signal,
    valid,
    baselines,
    thresholds,
):
    left = max(0, start - CONTEXT_SECONDS)
    right = min(len(signal), end + CONTEXT_SECONDS + 1)

    indices = np.arange(left, right)
    hours = indices / 3600

    event_start_hour = start / 3600
    event_end_hour = end / 3600

    fig, axes = plt.subplots(
        2, 1,
        figsize=(15, 10),
        sharex=True,
        gridspec_kw={"height_ratios": [1.3, 1]}
    )

    # --------------------------------------------------
    # TOP: RAW FLUX AND BASELINES
    # --------------------------------------------------

    axes[0].plot(
        hours,
        np.where(valid[left:right], signal[left:right], np.nan),
        label="Measured electron flux",
        color="black",
        linewidth=1.1,
    )

    for name, baseline in baselines.items():
        axes[0].plot(
            hours,
            baseline[left:right],
            label=name,
            linewidth=1.0,
        )

    axes[0].axvspan(
        event_start_hour,
        event_end_hour,
        alpha=0.15,
        color="grey",
        label="Reference event interval",
    )

    axes[0].set_ylabel(
        "Differential electron flux\n"
        r"(cm$^{-2}$ sr$^{-1}$ s$^{-1}$ keV$^{-1}$)"
    )

    axes[0].set_title(
        f"January 3, 2024 — Band 2 — Reference Event {event_number}\n"
        f"Samples {start}–{end} | "
        f"{event_start_hour:.3f}–{event_end_hour:.3f} hours"
    )

    axes[0].grid(True, alpha=0.3)
    axes[0].legend(loc="upper left", fontsize=8, ncol=2)

    # --------------------------------------------------
    # BOTTOM: RESIDUALS AND THRESHOLDS
    # --------------------------------------------------

    for name, baseline in baselines.items():
        residual = signal - baseline

        axes[1].plot(
            hours,
            residual[left:right],
            label=name,
            linewidth=1.0,
        )

        threshold = thresholds[name]["threshold"]

        axes[1].axhline(
            threshold,
            linestyle="--",
            linewidth=0.9,
            alpha=0.6,
        )

    axes[1].axvspan(
        event_start_hour,
        event_end_hour,
        alpha=0.15,
        color="grey",
    )

    axes[1].axhline(
        0,
        color="black",
        linewidth=0.8,
    )

    axes[1].set_xlabel("Time (hours from start of day)")
    axes[1].set_ylabel("Residual (flux units)")
    axes[1].set_title(
        "Residuals and method-specific positive detection thresholds"
    )

    axes[1].grid(True, alpha=0.3)
    axes[1].legend(loc="upper left", fontsize=8, ncol=2)

    fig.tight_layout()

    output = PLOT_DIR / (
        f"jan03_band2_event{event_number}_waveforms.png"
    )

    fig.savefig(output, dpi=180)
    plt.close(fig)

    print(f"Saved: {output}")


# --------------------------------------------------
# EVENT-LEVEL SUMMARY
# --------------------------------------------------

def generate_event_summary(
    signal,
    valid,
    baselines,
    thresholds,
):
    rows = []

    for event_number, (start, end) in enumerate(
        REFERENCE_EVENTS, start=1
    ):
        event_mask = np.zeros(len(signal), dtype=bool)
        event_mask[start:end + 1] = True

        event_mask &= valid

        for name, baseline in baselines.items():
            residual = signal - baseline

            threshold = thresholds[name]["threshold"]

            usable = event_mask & np.isfinite(residual)

            event_residual = residual[usable]
            active = event_residual > threshold

            rows.append({
                "event": event_number,
                "baseline": name,
                "threshold": threshold,
                "valid_samples": int(usable.sum()),
                "peak_residual": (
                    float(np.max(event_residual))
                    if len(event_residual) else np.nan
                ),
                "median_residual": (
                    float(np.median(event_residual))
                    if len(event_residual) else np.nan
                ),
                "active_samples": int(active.sum()),
                "active_percentage": (
                    100 * active.sum() / usable.sum()
                    if usable.sum() else np.nan
                ),
            })

    summary = pd.DataFrame(rows)

    output = (
        RESULT_DIR / "jan03_band2_event_waveform_summary.csv"
    )

    summary.to_csv(output, index=False)

    print("\nEVENT-LEVEL SUMMARY")
    print(summary.to_string(index=False))
    print(f"\nSaved summary: {output}")


# --------------------------------------------------
# MAIN
# --------------------------------------------------

def main():
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    signal, valid = load_data()

    baselines = calculate_baselines(signal)

    thresholds = calculate_thresholds(
        signal, valid, baselines
    )

    print("\nCALIBRATED THRESHOLDS")

    for name, values in thresholds.items():
        print(
            f"{name:28s} "
            f"threshold={values['threshold']:.2f} "
            f"robust_sigma={values['robust_sigma']:.2f}"
        )

    for event_number, (start, end) in enumerate(
        REFERENCE_EVENTS, start=1
    ):
        plot_event(
            event_number,
            start,
            end,
            signal,
            valid,
            baselines,
            thresholds,
        )

    generate_event_summary(
        signal,
        valid,
        baselines,
        thresholds,
    )

    print("\nAll waveform diagnostics completed.")


if __name__ == "__main__":
    main()
