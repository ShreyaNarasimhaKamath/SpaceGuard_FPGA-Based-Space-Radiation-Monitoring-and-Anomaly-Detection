
import numpy as np
from netCDF4 import Dataset

FILE = "data/raw/goes16_mpshi/goes16_mpshi_20240101.nc"

def decode_labels(variable):
    values = variable[:]
    labels = []

    for row in values:
        label = b"".join(row).decode("ascii", errors="replace")
        labels.append(label.strip().replace("\x00", ""))

    return labels


def print_stats(ds, name):
    var = ds.variables[name]
    values = np.ma.asarray(var[:], dtype=np.float64)
    values = np.ma.masked_invalid(values)

    valid = values.compressed()

    print(f"\n{name}")
    print("  Shape:", var.shape)
    print("  Units:", getattr(var, "units", "Not specified"))
    print("  Valid samples:", len(valid))

    if len(valid):
        print("  Minimum:", np.min(valid))
        print("  Maximum:", np.max(valid))
        print("  Mean:", np.mean(valid))
        print("  Median:", np.median(valid))

    for attr in ["_FillValue", "missing_value", "valid_min", "valid_max"]:
        if hasattr(var, attr):
            print(f"  {attr}:", getattr(var, attr))


with Dataset(FILE, "r") as ds:

    print("========== CHANNEL LABELS ==========")

    for name in [
        "direction_label",
        "energy_electron_label",
        "energy_proton_label"
    ]:
        print(f"\n{name}:")
        print(decode_labels(ds.variables[name]))

    print("\n========== SIGNAL STATISTICS ==========")

    for name in [
        "DiffElectronFluxes",
        "DiffProtonFluxes",
        "IntgElectronFluxes",
        "Dos1_HiLetDose",
        "Dos1_LoLetDose"
    ]:
        print_stats(ds, name)

    print("\n========== QUALITY FLAGS ==========")

    for name in [
        "DiffElectronFluxDQFs",
        "DiffProtonFluxDQFs",
        "IntgElectronFluxDQFs",
        "Dos1_HiLetDqf",
        "Dos1_LoLetDqf"
    ]:
        values, counts = np.unique(
            ds.variables[name][:],
            return_counts=True
        )
        print(f"\n{name}:")
        print(dict(zip(values.tolist(), counts.tolist())))

    print("\n========== TIMESTAMP CHECK ==========")

    timestamps = ds.variables["L1a_SciData_TimeStamp"][:]
    differences = np.diff(timestamps)

    print("First timestamp:", timestamps[0])
    print("Last timestamp:", timestamps[-1])
    print("Unique timestamp intervals:",
          np.unique(np.round(differences, 6), return_counts=True))

