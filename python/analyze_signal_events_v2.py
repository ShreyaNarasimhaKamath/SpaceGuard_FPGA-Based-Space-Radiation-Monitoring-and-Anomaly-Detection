
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset

FILE = "data/raw/goes16_mpshi/goes16_mpshi_20240101.nc"

# Configuration
BASELINE_WINDOW = 1800   # 30 minutes at 1 sample/second
CALIBRATION_HOURS = 10
THRESHOLD_MULTIPLIER = 6

with Dataset(FILE, "r") as ds:
    flux = np.asarray(
        ds.variables["DiffElectronFluxes"][:, 0, 0],
        dtype=np.float64
    )

    dqf = np.asarray(
        ds.variables["DiffElectronFluxDQFs"][:, 0, 0]
    )

    timestamps = np.asarray(
        ds.variables["L1a_SciData_TimeStamp"][:]
    )

# Keep only good-quality measurements
valid = (dqf == 0) & np.isfinite(flux) & (flux >= 0)

time_hours = (timestamps - timestamps[0]) / 3600

df = pd.DataFrame({
    "time": time_hours,
    "flux": np.where(valid, flux, np.nan)
})

# Causal rolling median: uses current and previous samples only
baseline = df["flux"].rolling(
    window=BASELINE_WINDOW,
    min_periods=BASELINE_WINDOW
).median()

residual = df["flux"] - baseline

# Initial 10-hour calibration interval
calibration = residual[
    (df["time"] < CALIBRATION_HOURS)
].dropna()

if len(calibration) == 0:
    raise RuntimeError("No valid calibration residuals.")

# Robust noise estimate
residual_median = calibration.median()

mad = np.median(
    np.abs(calibration.to_numpy() - residual_median)
)

robust_sigma = 1.4826 * mad

threshold = THRESHOLD_MULTIPLIER * robust_sigma

# Fixed threshold; does not increase during candidate events
candidate = residual > threshold

print("\n========== DETECTOR V2 ==========")
print("Baseline window:", BASELINE_WINDOW, "seconds")
print("Calibration interval:", CALIBRATION_HOURS, "hours")
print("Calibration samples:", len(calibration))
print("Calibration residual median:", residual_median)
print("Robust sigma:", robust_sigma)
print("Detection threshold:", threshold)
print("Candidate samples:", int(candidate.sum()))

# Group consecutive candidate samples
active = candidate.fillna(False).to_numpy()
times = df["time"].to_numpy()

starts = np.where(active & ~np.r_[False, active[:-1]])[0]
ends = np.where(active & ~np.r_[active[1:], False])[0]

print("\nCandidate periods:")

for start, end in zip(starts, ends):
    print(
        f"{times[start]:.3f} h to {times[end]:.3f} h | "
        f"{end - start + 1} seconds"
    )

# Plot signal and causal baseline
plt.figure(figsize=(14, 6))
plt.plot(df["time"], df["flux"], linewidth=0.5,
         label="Measured flux")
plt.plot(df["time"], baseline, linewidth=1,
         label="Causal 30-minute median")

# Group nearby threshold crossings into candidate events

# Group nearby threshold crossings into candidate events

active = candidate.fillna(False).to_numpy()
times = df["time"].to_numpy()
flux_values = df["flux"].to_numpy()

MAX_GAP_SECONDS = 60
MIN_ACTIVE_SECONDS = 60

indices = np.where(active)[0]
events = []

if len(indices):
    start = previous = indices[0]
    active_count = 1

    for idx in indices[1:]:
        gap = idx - previous - 1

        if gap <= MAX_GAP_SECONDS:
            active_count += 1
            previous = idx
        else:
            if active_count >= MIN_ACTIVE_SECONDS:
                events.append((start, previous, active_count))

            start = previous = idx
            active_count = 1

    if active_count >= MIN_ACTIVE_SECONDS:
        events.append((start, previous, active_count))

print("\n========== GROUPED CANDIDATE EVENTS ==========")
print("Maximum merge gap:", MAX_GAP_SECONDS, "seconds")
print("Minimum active duration:", MIN_ACTIVE_SECONDS, "seconds")
print("Number of candidate events:", len(events))

for number, (start, end, active_count) in enumerate(events, 1):
    span = end - start + 1
    peak = np.nanmax(flux_values[start:end + 1])

    print(
        f"Event {number}: "
        f"{times[start]:.3f} h to {times[end]:.3f} h | "
        f"Span: {span} seconds | "
        f"Active samples: {active_count} | "
        f"Peak flux: {peak:.2f}"
    )




# Merge detections separated by gaps of at most 60 seconds
MAX_GAP_SECONDS = 60

# Reject events shorter than 60 seconds
MIN_EVENT_SECONDS = 60
indices = np.where(active)[0]
events = []

if len(indices):
    start = previous = indices[0]

    for idx in indices[1:]:
        gap = idx - previous - 1

        if gap <= MAX_GAP_SECONDS:
            previous = idx
        else:
            duration = previous - start + 1

            if duration >= MIN_EVENT_SECONDS:
                events.append((start, previous))

            start = previous = idx

    duration = previous - start + 1

    if duration >= MIN_EVENT_SECONDS:
        events.append((start, previous))

print("\n========== GROUPED CANDIDATE EVENTS ==========")
print("Maximum merge gap:", MAX_GAP_SECONDS, "seconds")
print("Minimum event duration:", MIN_EVENT_SECONDS, "seconds")
print("Number of candidate events:", len(events))

for number, (start, end) in enumerate(events, 1):
    duration = end - start + 1
    peak = df["flux"].iloc[start:end + 1].max()

    print(
        f"Event {number}: "
        f"{times[start]:.3f} h to {times[end]:.3f} h | "
        f"Duration: {duration} seconds | "
        f"Peak flux: {peak:.2f}"
    )

plt.xlabel("Time (hours)")
plt.ylabel("Differential electron flux")
plt.title("Causal Baseline Estimation")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("results/plots/causal_baseline.png", dpi=200)

# Plot residual and fixed threshold
plt.figure(figsize=(14, 6))
plt.plot(df["time"], residual, linewidth=0.5,
         label="Flux residual")
plt.axhline(threshold, linestyle="--",
            label="Fixed positive threshold")

plt.xlabel("Time (hours)")
plt.ylabel("Flux deviation")
plt.title("SpaceGuard Detector V2")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("results/plots/anomaly_candidates_v2.png", dpi=200)

plt.show()
