
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

FILE = "data/raw/goes16_mpshi/goes16_mpshi_20240101.nc"

START = 11.475
END = 13.25

with Dataset(FILE, "r") as ds:
    flux = np.asarray(
        ds.variables["DiffElectronFluxes"][:],
        dtype=np.float64
    )
    dqf = np.asarray(
        ds.variables["DiffElectronFluxDQFs"][:]
    )
    timestamps = np.asarray(
        ds.variables["L1a_SciData_TimeStamp"][:]
    )

time = (timestamps - timestamps[0]) / 3600
mask = (time >= START) & (time <= END)

fig, axes = plt.subplots(
    3, 1, figsize=(14, 11), sharex=True
)

for band in range(3):
    ax = axes[band]

    for telescope in range(5):
        values = flux[mask, telescope, band].copy()
        flags = dqf[mask, telescope, band]

        values[flags != 0] = np.nan
        values[~np.isfinite(values)] = np.nan
        values[values < 0] = np.nan

        ax.plot(
            time[mask],
            values,
            linewidth=0.8,
            label=f"Telescope {telescope + 1}"
        )

    ax.set_title(f"Electron Energy Band {band + 1}")
    ax.set_ylabel("Differential flux")
    ax.grid(True)
    ax.legend(ncol=5)

axes[-1].set_xlabel("Time (hours)")
fig.suptitle(
    "Event 2: Energy-Band Comparison Across Five Telescopes"
)

plt.tight_layout()
plt.savefig(
    "results/plots/event2_energy_comparison.png",
    dpi=200
)
plt.show()
