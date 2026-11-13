
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date

FILE = "data/raw/goes16_mpshi/goes16_mpshi_20240101.nc"

with Dataset(FILE, "r") as ds:

    flux = ds.variables["DiffElectronFluxes"][:, 0, 0]
    dqf = ds.variables["DiffElectronFluxDQFs"][:, 0, 0]
    timestamps = ds.variables["L1a_SciData_TimeStamp"][:]

    flux = np.asarray(flux, dtype=np.float64)
    dqf = np.asarray(dqf)

    # Strict clean-data selection
    valid = (dqf == 0) & np.isfinite(flux) & (flux >= 0)

    clean_flux = flux[valid]

    print("Selected channel: Telescope 1, Electron Band 1")
    print("Total samples:", len(flux))
    print("Clean samples:", len(clean_flux))
    print("Rejected samples:", len(flux) - len(clean_flux))

    print("\nSignal statistics:")
    for p in [0, 25, 50, 75, 90, 95, 99, 99.9, 100]:
        print(f"{p:5.1f} percentile: {np.percentile(clean_flux, p):.4f}")

    # Convert timestamps to hours from beginning of day
    time_hours = (timestamps - timestamps[0]) / 3600

    # Plot full-day signal
    plt.figure(figsize=(14, 5))
    plt.plot(time_hours[valid], flux[valid], linewidth=0.6)
    plt.xlabel("Time (hours from start of day)")
    plt.ylabel("Differential Electron Flux")
    plt.title("GOES-16 MPS-HI: Telescope 1, Electron Band 1")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("results/plots/electron_band1_full_day.png", dpi=200)

    # Plot logarithmic view
    plt.figure(figsize=(14, 5))
    plt.plot(time_hours[valid], clean_flux, linewidth=0.6)
    plt.yscale("log")
    plt.xlabel("Time (hours from start of day)")
    plt.ylabel("Differential Electron Flux (log scale)")
    plt.title("Electron Band 1 - Logarithmic View")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("results/plots/electron_band1_log.png", dpi=200)

    # Plot first hour
    first_hour = time_hours < 1

    plt.figure(figsize=(14, 5))
    plt.plot(time_hours[first_hour & valid],
             flux[first_hour & valid], linewidth=0.8)
    plt.xlabel("Time (hours from start of day)")
    plt.ylabel("Differential Electron Flux")
    plt.title("Electron Band 1 - First Hour")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("results/plots/electron_band1_first_hour.png", dpi=200)

    plt.show()

print("\nPlots saved in results/plots/")

