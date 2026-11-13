
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

FILE = "data/raw/goes16_mpshi/goes16_mpshi_20240101.nc"

EVENT_START = 11.975
EVENT_END = 12.750

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

time_hours = (timestamps - timestamps[0]) / 3600

mask = (
    (time_hours >= EVENT_START - 0.5)
    & (time_hours <= EVENT_END + 0.5)
)

fig, axes = plt.subplots(5, 1, figsize=(14, 13), sharex=True)

for telescope in range(5):

    ax = axes[telescope]

    for energy in range(3):

        values = flux[mask, telescope, energy].copy()
        flags = dqf[mask, telescope, energy]

        values[flags != 0] = np.nan

        ax.plot(
            time_hours[mask],
            values,
            linewidth=0.8,
            label=f"Electron Band {energy + 1}"
        )

    ax.set_ylabel("Flux")
    ax.set_title(f"Telescope {telescope + 1}")
    ax.legend(loc="upper right")
    ax.grid(True)

axes[-1].set_xlabel("Time (hours)")
fig.suptitle(
    "SpaceGuard: Cross-Telescope Event 2 Analysis",
    fontsize=14
)

plt.tight_layout()
plt.savefig(
    "results/plots/event2_multichannel.png",
    dpi=200
)
plt.show()
