
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from pathlib import Path

# --------------------------------------------------
# CONFIGURATION
# --------------------------------------------------

FILE = "data/raw/goes16_mpshi/goes16_mpshi_20240101.nc"

TELESCOPE = 0
ENERGY_BANDS = [0, 1]       # Electron Bands 1 and 2

BASELINE_WINDOW = 1800     # 30 minutes
CALIBRATION_HOURS = 10
THRESHOLD_MULTIPLIER = 6

MAX_GAP_SECONDS = 60       # Merge nearby detections
MIN_EVENT_SECONDS = 60     # Minimum event duration
MATCH_TOLERANCE_SECONDS = 300  # 5-minute matching tolerance

OUTPUT = Path("results/evaluation")
PLOTS = Path("results/plots")

OUTPUT.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------
# EVENT GROUPING
# --------------------------------------------------

def group_events(candidate, times, flux):
    active = np.asarray(candidate, dtype=bool)
    indices = np.where(active)[0]

    events = []

    if len(indices) == 0:
        return events

    start = previous = indices[0]

    groups = []

    for idx in indices[1:]:
        gap = idx - previous - 1

        if gap <= MAX_GAP_SECONDS:
            previous = idx
        else:
            groups.append((start, previous))
            start = previous = idx

    groups.append((start, previous))

    for start, end in groups:
        duration = end - start + 1

        if duration < MIN_EVENT_SECONDS:
            continue

        events.append({
            "start_index": int(start),
            "end_index": int(end),
            "start_hour": float(times[start]),
            "end_hour": float(times[end]),
            "duration_seconds": int(duration),
            "peak_flux": float(np.nanmax(flux[start:end + 1]))
        })

    return events


# --------------------------------------------------
# LOAD DATA
# --------------------------------------------------

with Dataset(FILE, "r") as ds:

    timestamps = np.asarray(
        ds.variables["L1a_SciData_TimeStamp"][:],
        dtype=np.float64
    )

    channel_data = {}

    for band in ENERGY_BANDS:

        flux = np.asarray(
            ds.variables["DiffElectronFluxes"][:, TELESCOPE, band],
            dtype=np.float64
        )

        dqf = np.asarray(
            ds.variables["DiffElectronFluxDQFs"][:, TELESCOPE, band]
        )

        valid = (
            (dqf == 0)
            & np.isfinite(flux)
            & (flux >= 0)
        )

        channel_data[band] = {
            "flux": flux,
            "valid": valid
        }


time_hours = (timestamps - timestamps[0]) / 3600

if len(timestamps) < BASELINE_WINDOW:
    raise RuntimeError("Dataset is shorter than the baseline window.")


# --------------------------------------------------
# INDEPENDENT CHANNEL PROCESSING
# --------------------------------------------------

results = {}

for band in ENERGY_BANDS:

    flux = channel_data[band]["flux"]
    valid = channel_data[band]["valid"]

    signal = pd.Series(
        np.where(valid, flux, np.nan)
    )

    # Causal rolling baseline
    baseline = signal.rolling(
        window=BASELINE_WINDOW,
        min_periods=BASELINE_WINDOW
    ).median()

    residual = signal - baseline

    # Initial calibration interval
    calibration = residual[
        (time_hours < CALIBRATION_HOURS)
    ].dropna()

    if len(calibration) == 0:
        raise RuntimeError(
            f"No calibration samples for Band {band + 1}"
        )

    calibration_values = calibration.to_numpy()

    residual_median = np.median(calibration_values)

    mad = np.median(
        np.abs(calibration_values - residual_median)
    )

    robust_sigma = 1.4826 * mad

    if not np.isfinite(robust_sigma) or robust_sigma <= 0:
        raise RuntimeError(
            f"Invalid robust noise estimate for Band {band + 1}"
        )

    threshold = THRESHOLD_MULTIPLIER * robust_sigma

    # Positive deviations only
    candidate = (
        (residual > threshold)
        & valid
        & baseline.notna()
    ).to_numpy()

    events = group_events(
        candidate,
        time_hours,
        flux
    )

    results[band] = {
        "flux": flux,
        "valid": valid,
        "baseline": baseline.to_numpy(),
        "residual": residual.to_numpy(),
        "threshold": threshold,
        "robust_sigma": robust_sigma,
        "calibration_samples": len(calibration),
        "candidate_samples": int(candidate.sum()),
        "events": events
    }

    print(f"\n========== ELECTRON BAND {band + 1} ==========")
    print("Valid samples:", int(valid.sum()))
    print("Calibration samples:", len(calibration))
    print("Robust sigma:", robust_sigma)
    print("Detection threshold:", threshold)
    print("Candidate samples:", int(candidate.sum()))
    print("Detected events:", len(events))

    for i, event in enumerate(events, 1):
        print(
            f"Event {i}: "
            f"{event['start_hour']:.3f} - "
            f"{event['end_hour']:.3f} hours | "
            f"Duration: {event['duration_seconds']} s | "
            f"Peak: {event['peak_flux']:.2f}"
        )


# --------------------------------------------------
# CROSS-CHANNEL EVENT MATCHING
# --------------------------------------------------

band1_events = results[0]["events"]
band2_events = results[1]["events"]

tolerance_hours = MATCH_TOLERANCE_SECONDS / 3600

matched_pairs = []

for i, event1 in enumerate(band1_events):

    for j, event2 in enumerate(band2_events):

        # Distance between two event intervals
        gap = max(
            0,
            event1["start_hour"] - event2["end_hour"],
            event2["start_hour"] - event1["end_hour"]
        )

        if gap <= tolerance_hours:
            matched_pairs.append((i, j))

matched_band1 = {i for i, j in matched_pairs}
matched_band2 = {j for i, j in matched_pairs}

print("\n========== CROSS-CHANNEL MATCHING ==========")
print("Matching tolerance:", MATCH_TOLERANCE_SECONDS, "seconds")
print("Matched event pairs:", len(matched_pairs))
print("Band 1-only events:", len(band1_events) - len(matched_band1))
print("Band 2-only events:", len(band2_events) - len(matched_band2))

for i, j in matched_pairs:
    print(
        f"Band 1 event {i + 1} <-> Band 2 event {j + 1}"
    )


# --------------------------------------------------
# SAVE EVENT RESULTS
# --------------------------------------------------

rows = []

for band, channel in results.items():

    for number, event in enumerate(channel["events"], 1):

        rows.append({
            "energy_band": band + 1,
            "event_number": number,
            **event
        })

event_df = pd.DataFrame(rows)

event_df.to_csv(
    OUTPUT / "independent_channel_events.csv",
    index=False
)

print("\nEvent results saved to:")
print(OUTPUT / "independent_channel_events.csv")


# --------------------------------------------------
# PLOT SIGNALS AND BASELINES
# --------------------------------------------------

fig, axes = plt.subplots(
    2, 1,
    figsize=(15, 9),
    sharex=True
)

for ax, band in zip(axes, ENERGY_BANDS):

    channel = results[band]

    ax.plot(
        time_hours,
        channel["flux"],
        linewidth=0.45,
        label="Measured flux"
    )

    ax.plot(
        time_hours,
        channel["baseline"],
        linewidth=1,
        label="Causal baseline"
    )

    ax.set_title(f"Electron Band {band + 1}")
    ax.set_ylabel("Differential flux")
    ax.grid(True)
    ax.legend()

axes[-1].set_xlabel("Time (hours)")

plt.tight_layout()
plt.savefig(
    PLOTS / "independent_channel_baselines.png",
    dpi=200
)


# --------------------------------------------------
# PLOT RESIDUALS AND THRESHOLDS
# --------------------------------------------------

fig, axes = plt.subplots(
    2, 1,
    figsize=(15, 9),
    sharex=True
)

for ax, band in zip(axes, ENERGY_BANDS):

    channel = results[band]

    ax.plot(
        time_hours,
        channel["residual"],
        linewidth=0.45,
        label="Flux residual"
    )

    ax.axhline(
        channel["threshold"],
        linestyle="--",
        label="Positive detection threshold"
    )

    ax.set_title(f"Electron Band {band + 1} Detection")
    ax.set_ylabel("Flux deviation")
    ax.grid(True)
    ax.legend()

axes[-1].set_xlabel("Time (hours)")

plt.tight_layout()
plt.savefig(
    PLOTS / "independent_channel_residuals.png",
    dpi=200
)

plt.show()
