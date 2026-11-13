from netCDF4 import Dataset

file_path = "data/raw/goes16_mpshi/goes16_mpshi_20240101.nc"

with Dataset(file_path, "r") as ds:

    print("\n========== DATASET INFORMATION ==========")

    print("Dimensions:")
    for name, dimension in ds.dimensions.items():
        print(name, len(dimension))

    print("\nVariables:")
    for name, variable in ds.variables.items():
        print("\nVariable:", name)
        print("Dimensions:", variable.dimensions)
        print("Shape:", variable.shape)
        print("Data type:", variable.dtype)
        print("Units:", getattr(variable, "units", "Not specified"))
        print("Description:", getattr(variable, "long_name", "Not specified"))

    print("\nGlobal attributes:")
    for name in ds.ncattrs():
        print(name, ":", getattr(ds, name))
