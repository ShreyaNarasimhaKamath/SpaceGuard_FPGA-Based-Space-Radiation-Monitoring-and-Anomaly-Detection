
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import pandas as pd

FILE = "data/raw/goes16_mpshi/goes16_mpshi_20240101.nc"

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

# Retain only good-quality measurements
valid = (dqf == 0) & np.isfinite(flux) & (flux >= 0)

df = pd.DataFrame({
    "time": (timestamps - timestamps[0]) / 3600,
    "flux": flux
})

df.loc[~valid, "flux"] = np.nan

# 5-minute rolling baseline
baseline = df["flux"].rolling(
    window=300,
    min_periods=300,
    center=True
).median()

# Deviation from baseline
deviation = df["flux"] - baseline

# Robust rolling variability estimate
absolute_deviation = (df["flux"] - baseline).abs()

rolling_mad = absolute_deviation.rolling(
    window=300,
    min_periods=300,
    center=True
).median()

# Candidate threshold: 6 times rolling MAD
threshold = 6 * 1.4826 * rolling_mad

candidate = deviation > threshold

print("Total samples:", len(df))
print("Valid samples:", valid.sum())
print("Candidate anomaly samples:", candidate.sum())

print("\nCandidate periods:")

active = candidate.fillna(False).to_numpy()
times = df["time"].to_numpy()

starts = np.where(active & ~np.r_[False, active[:-1]])[0]
ends = np.where(active & ~np.r_[active[1:], False])[0]

for start, end in zip(starts, ends):
    print(
        f"{times[start]:.3f} h to "
        f"{times[end]:.3f} h | "
        f"{end - start + 1} seconds"
    )

plt.figure(figsize=(14, 6))
plt.plot(df["time"], df["flux"], label="Measured flux", linewidth=0.6)
plt.plot(df["time"], baseline, label="5-minute rolling median")
plt.xlabel("Time (hours)")
plt.ylabel("Differential electron flux")
plt.title("Radiation Signal and Rolling Baseline")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("results/plots/rolling_baseline.png", dpi=200)

plt.figure(figsize=(14, 5))
plt.plot(df["time"], deviation, label="Baseline deviation")
plt.plot(df["time"], threshold, label="Positive detection threshold")
plt.xlabel("Time (hours)")
plt.ylabel("Flux deviation")
plt.title("Candidate Anomaly Detection")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("results/plots/anomaly_candidates.png", dpi=200)

plt.show()
