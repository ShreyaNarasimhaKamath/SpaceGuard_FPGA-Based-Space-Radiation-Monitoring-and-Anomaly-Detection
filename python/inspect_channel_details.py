
import numpy as np
from netCDF4 import Dataset

FILE = "data/raw/goes16_mpshi/goes16_mpshi_20240101.nc"

with Dataset(FILE, "r") as ds:

    print("========== DQF VARIABLE ATTRIBUTES ==========")

    for name in [
        "DiffElectronFluxDQFs",
        "DiffProtonFluxDQFs",
        "Dos1_HiLetDqf",
        "Dos1_LoLetDqf"
    ]:
        var = ds.variables[name]

        print(f"\n{name}")
        for attr in var.ncattrs():
            print(f"  {attr}: {getattr(var, attr)}")

    print("\n========== ELECTRON CHANNEL STATISTICS ==========")

    flux = ds.variables["DiffElectronFluxes"]
    dqf = ds.variables["DiffElectronFluxDQFs"]

    for direction in range(5):
        for energy in range(10):

            x = np.ma.asarray(
                flux[:, direction, energy],
                dtype=np.float64
            )

            q = np.asarray(
                dqf[:, direction, energy]
            )

            valid = x.compressed()

            if len(valid) == 0:
                continue

            print(
                f"Telescope {direction + 1}, "
                f"Electron {energy + 1}: "
                f"min={valid.min():.5g}, "
                f"max={valid.max():.5g}, "
                f"mean={valid.mean():.5g}, "
                f"negative={(valid < 0).sum()}, "
                f"DQF={dict(zip(*np.unique(q, return_counts=True)))}"
            )

    print("\n========== PROTON CHANNEL STATISTICS ==========")

    flux = ds.variables["DiffProtonFluxes"]
    dqf = ds.variables["DiffProtonFluxDQFs"]

    for direction in range(5):
        for energy in range(11):

            x = np.ma.asarray(
                flux[:, direction, energy],
                dtype=np.float64
            )

            q = np.asarray(
                dqf[:, direction, energy]
            )

            valid = x.compressed()

            if len(valid) == 0:
                continue

            print(
                f"Telescope {direction + 1}, "
                f"Proton {energy + 1}: "
                f"min={valid.min():.5g}, "
                f"max={valid.max():.5g}, "
                f"mean={valid.mean():.5g}, "
                f"negative={(valid < 0).sum()}, "
                f"DQF={dict(zip(*np.unique(q, return_counts=True)))}"
            )

