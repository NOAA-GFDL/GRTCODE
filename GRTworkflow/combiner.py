import os
import glob
import xarray as xr
import numpy as np
import argparse

def combine_netcdf_files(input_files, output_file):
    """
    Combine NetCDF files by integrating selected variables over 'lw_wavenumber'
    and concatenating along the 'lon' dimension. Only variables starting with 'rl'
    and listed in 'variables_to_process' are included.
    """
    variables_to_process = ['rlutcsaf', 'rldscsaf', 'rlutaf', 'rlucsaf_level', 'rluaf_level', 'rldcsaf_level', 'rldaf_level']
    datasets = []

    for file in input_files:
        print(f"Processing file: {file}")
        ds = xr.open_dataset(file)

        integrated_vars = {}

        for var in variables_to_process:
            if var in ds.data_vars:
                data = ds[var]

                # Replace NaN _FillValue if needed
                if '_FillValue' in data.attrs and np.isnan(data.attrs['_FillValue']):
                    data.attrs['_FillValue'] = -9999.0

                # Convert to float32 if needed
                if data.dtype == 'float64':
                    data = data.astype('float32')
                # Rebin from 0.1 cm-1 to 1 cm-1
                if 'lw_wavenumber' in data.dims:

                    # sanity check (optional but recommended)
                    dnu = float(data['lw_wavenumber'].diff('lw_wavenumber').mean())
                    if not np.isclose(dnu, 0.1):
                        raise ValueError(f"Unexpected spectral resolution: {dnu}")

                    rebinned = (
                        data
                        .coarsen(lw_wavenumber=10, boundary='trim')
                        .sum()
                    )

                    # redefine spectral coordinate as bin centers
                    new_nu = (
                        data['lw_wavenumber']
                        .coarsen(lw_wavenumber=10, boundary='trim')
                        .mean()
                    )

                    rebinned = rebinned.assign_coords(lw_wavenumber=new_nu)

                    rebinned.name = f"{var}"
                    rebinned.attrs['units'] = "W m-2"
                    rebinned.attrs['long_name'] = f"{var} rebinned to 1 cm-1"

                    integrated_vars[rebinned.name] = rebinned / 10.0
                else:
                    print(f"Skipping {var} — missing 'lw_wavenumber' dimension")
            else:
                print(f"Variable '{var}' not found in: {file}")
        # Create dataset with integrated variables
        coords_to_keep = {dim: ds[dim] for dim in ds.dims if dim != 'lw_wavenumber'}
        integrated_ds = xr.Dataset(integrated_vars, coords=coords_to_keep)
        datasets.append(integrated_ds)

    # Concatenate along 'lon' and sort
    combined_ds = xr.concat(
        datasets,
        dim='lon',
        data_vars='all',
        coords='all',
        compat='override',
        combine_attrs='drop'
    ).sortby('lon')

    # Add longitude metadata
    combined_ds['lon'].attrs.update({
        'long_name': 'longitude',
        'units': 'degrees_E',
        'axis': 'X'
    })

    # Save to NetCDF
    combined_ds.to_netcdf(output_file)
    print("Final longitude ordering:", combined_ds['lon'].values)
    print(f"Integrated NetCDF file saved to {output_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--workdir', type=str, required=True, help='Path to the working directory')
    parser.add_argument('--year', type=int, required=True, help='Year of the data to combine')
    parser.add_argument('--nametag', type=str, required=True, help='Experiment name')
    args = parser.parse_args()

    file_pattern = f"{args.workdir}/{args.nametag}/{args.year}.era5-fluxes.nc.*"
    output_file = f"{args.workdir}/{args.year}.{args.nametag}-cleansky-integrated.nc"

    input_files = sorted(glob.glob(file_pattern))
    if not input_files:
        print(f"No files found matching pattern: {file_pattern}")
    else:
        combine_netcdf_files(input_files, output_file)
