
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset

FILE = "data/raw/goes16_mpshi/goes16_mpshi_20240101.nc"

# Detector configuration
BASELINE_WINDOW = 1800
CALIBRATION_HOURS = 10
THRESHOLD_MULTIPLIER = 6

MAX_GAP_SECONDS = 60
MIN_ACTIVE_SECONDS = 60

with Dataset(FILE, "r") as ds:

    electron = np.asarray(
        ds.variables["DiffElectronFluxes"][:],
        dtype=np.float64
    )

    electron_dqf = np.asarray(
        ds.variables["DiffElectronFluxDQFs"][:]
    )

    timestamps = np.asarray(
        ds.variables["L1a_SciData_TimeStamp"][:]
    )

time_hours = (timestamps - timestamps[0]) / 3600

# --------------------------------------------------
# 1. Reconstruct the reference detector
# --------------------------------------------------

reference_flux = electron[:, 0, 0]
reference_dqf = electron_dqf[:, 0, 0]

valid = (
    (reference_dqf == 0)
    & np.isfinite(reference_flux)
    & (reference_flux >= 0)
)

reference = pd.Series(
    np.where(valid, reference_flux, np.nan)
)

baseline = reference.rolling(
    window=BASELINE_WINDOW,
    min_periods=BASELINE_WINDOW
).median()

residual = reference - baseline

calibration = residual[
    np.asarray(time_hours < CALIBRATION_HOURS)
].dropna()

residual_median = calibration.median()

mad = np.median(
    np.abs(calibration.to_numpy() - residual_median)
)

robust_sigma = 1.4826 * mad
threshold = THRESHOLD_MULTIPLIER * robust_sigma

candidate = (residual > threshold).fillna(False).to_numpy()

# --------------------------------------------------
# 2. Group candidate samples into events
# --------------------------------------------------

indices = np.where(candidate)[0]
events = []

if len(indices):
    start = previous = indices[0]
    active_count = 1

    for idx in indices[1:]:
        gap = idx - previous - 1

        if gap <= MAX_GAP_SECONDS:
            previous = idx
            active_count += 1
        else:
            if active_count >= MIN_ACTIVE_SECONDS:
                events.append((start, previous))

            start = previous = idx
            active_count = 1

    if active_count >= MIN_ACTIVE_SECONDS:
        events.append((start, previous))

print("\n========== REFERENCE DETECTOR ==========")
print("Threshold:", threshold)
print("Candidate events:", len(events))

# --------------------------------------------------
# 3. Cross-channel event validation
# --------------------------------------------------

results = []

for event_number, (start, end) in enumerate(events, 1):

    # Use the preceding 30 minutes as comparison interval
    pre_start = max(0, start - BASELINE_WINDOW)
    pre_end = start

    event_flux = electron[start:end + 1]
    pre_flux = electron[pre_start:pre_end]

    event_dqf = electron_dqf[start:end + 1]
    pre_dqf = electron_dqf[pre_start:pre_end]

    for direction in range(5):
        for energy in range(10):

            event_values = event_flux[:, direction, energy]
            event_flags = event_dqf[:, direction, energy]

            pre_values = pre_flux[:, direction, energy]
            pre_flags = pre_dqf[:, direction, energy]

            event_valid = (
                (event_flags == 0)
                & np.isfinite(event_values)
                & (event_values >= 0)
            )

            pre_valid = (
                (pre_flags == 0)
                & np.isfinite(pre_values)
                & (pre_values >= 0)
            )

            ev = event_values[event_valid]
            pv = pre_values[pre_valid]

            if len(ev) < 10 or len(pv) < 10:
                continue

            event_median = np.median(ev)
            pre_median = np.median(pv)

            if pre_median <= 0:
                continue

            fold_change = event_median / pre_median
            percent_change = (fold_change - 1) * 100

            results.append({
                "event": event_number,
                "start_hour": time_hours[start],
                "end_hour": time_hours[end],
                "telescope": direction + 1,
                "energy_band": energy + 1,
                "event_median": event_median,
                "pre_event_median": pre_median,
                "fold_change": fold_change,
                "percent_change": percent_change,
                "event_valid_samples": len(ev),
                "event_total_samples": len(event_values)
            })

# --------------------------------------------------
# 4. Display results
# --------------------------------------------------

df_results = pd.DataFrame(results)

if df_results.empty:
    raise RuntimeError("No valid cross-channel comparisons found.")

print("\n========== CROSS-CHANNEL RESULTS ==========")

for event_number in sorted(df_results["event"].unique()):

    subset = df_results[
        df_results["event"] == event_number
    ]

    print(f"\nEVENT {event_number}")
    print("Time:",
          f"{subset['start_hour'].iloc[0]:.3f} - "
          f"{subset['end_hour'].iloc[0]:.3f} hours")

    print("Channels showing the largest increases:")

    top = subset.sort_values(
        "fold_change", ascending=False
    ).head(8)

    print(
        top[
            [
                "telescope",
                "energy_band",
                "fold_change",
                "percent_change"
            ]
        ].to_string(index=False)
    )

# Save all comparisons
output = "results/evaluation/cross_channel_results.csv"

df_results.to_csv(output, index=False)

print("\nResults saved to:", output)

# --------------------------------------------------
# 5. Summary plot
# --------------------------------------------------

summary = df_results.groupby("event")["percent_change"].median()

plt.figure(figsize=(12, 5))
summary.plot(kind="bar")

plt.xlabel("Candidate event")
plt.ylabel("Median cross-channel change (%)")
plt.title("Cross-Channel Particle Flux Variation")
plt.grid(axis="y")
plt.tight_layout()

plt.savefig(
    "results/plots/cross_channel_summary.png",
    dpi=200
)

plt.show()
