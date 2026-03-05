import os
import numpy as np
from netCDF4 import Dataset
import subprocess

# Path containing the files
path = "./your_directory_path"  # Replace with your directory path
output_path = "./broadband"  # Folder to save the output files
os.makedirs(output_path, exist_ok=True)  # Create the folder if it doesn't exist

# Define the experiment name
exp_name = "PI"

# Define the range of years
years = range(2003, 2024)

# Iterate through each year
for year in years:
    file_name = f"{year}.{exp_name}.nc"
    file_path = os.path.join(path, file_name)

    # Check if the file exists
    if not os.path.exists(file_path):
        print(f"File does not exist: {file_name}")
        continue

    # Open the NetCDF file and process the variables
    try:
        with Dataset(file_path, "r") as nc_file:
            variables_to_process = ['rlutcsaf', 'rldscsaf', 'rlutaf', 'rlutropcsaf', 'rlutropaf', 'rldtropcsaf', 'rldtropaf']
            integrated_data = {}

            for var_name in variables_to_process:
                if var_name in nc_file.variables:
                    var_data = nc_file.variables[var_name][:]
                    integrated_data[var_name] = var_data.sum(axis=-1)
                else:
                    print(f"Variable '{var_name}' not found in: {file_name}")

            # Save integrated values to a new NetCDF file
            output_file_path = os.path.join(output_path, file_name)
            with Dataset(output_file_path, "w") as output_nc:
                # Copy dimensions excluding wavenumber
                for dim_name, dim in nc_file.dimensions.items():
                    if dim_name != 'wavenumber':
                        output_nc.createDimension(dim_name, len(dim) if not dim.isunlimited() else None)

                # Copy latitude and longitude variables if they exist
                for var_name in ['latitude', 'longitude']:
                    if var_name in nc_file.variables:
                        var = nc_file.variables[var_name]
                        new_var = output_nc.createVariable(var_name, var.datatype, var.dimensions)
                        new_var[:] = var[:]
                        new_var.setncatts({attr: var.getncattr(attr) for attr in var.ncattrs()})

                # Create new variables for the integrated values
                for var_name, integrated_var in integrated_data.items():
                    integrated_var_name = f"{var_name}_integrated"
                    integrated_var_nc = output_nc.createVariable(integrated_var_name, integrated_var.dtype, ('time', 'lat', 'lon'))
                    integrated_var_nc[:] = integrated_var
                    integrated_var_nc.units = "Wm^-2"
                    integrated_var_nc.long_name = f"Spectrally integrated {var_name}"

            print(f"Saved summed variables to {output_file_path}")
    except Exception as e:
        print(f"Error processing {file_name}: {e}")

