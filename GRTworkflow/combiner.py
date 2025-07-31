import os
import glob
import xarray as xr
import argparse

def combine_netcdf_files(input_files, output_file):
    """
    Combine NetCDF files by concatenating along the 'lon' dimension.
    Only variables starting with 'rl' are included, and longitude is globally sorted.
    """
    datasets = []
    for file in input_files:
        ds = xr.open_dataset(file)

        # Filter variables starting with 'rl'
        rl_vars = {}
        for var in ds.data_vars:
            if var.startswith('rl'):
                data = ds[var]
                if '_FillValue' in data.attrs and np.isnan(data.attrs['_FillValue']):
                    data.attrs['_FillValue'] = -9999.0
                if data.dtype == 'float64':
                    data = data.astype('float32')
                rl_vars[var] = data

        rl_ds = xr.Dataset(rl_vars, coords=ds.coords)
        datasets.append(rl_ds)

    # Concatenate and globally sort by 'lon'
    combined_ds = xr.concat(
        datasets,
        dim='lon',
        data_vars='all',
        coords='all',
        compat='override',
        combine_attrs='drop'
    ).sortby('lon')

    # Save to NetCDF
    combined_ds['lon'].attrs['long_name'] = 'longitude'
    combined_ds['lon'].attrs['units'] = 'degrees_E'
    combined_ds['lon'].attrs['axis'] = 'X'

    # --- Dimension reordering section ---
    target_order = ("time", "lw_wavenumber", "lat", "lon")
    reordered_vars = {}

    for varname in combined_ds.data_vars:
        var = combined_ds[varname]

        if var.ndim <= 1:
            reordered_vars[varname] = var
            continue

        if all(dim in var.dims for dim in target_order):
            reordered_vars[varname] = var.transpose(*target_order)
        else:
            print(f"Skipping variable '{varname}' — unexpected dimensions: {var.dims}")
            reordered_vars[varname] = var  # Keep as-is

    final_ds = xr.Dataset(reordered_vars, coords=combined_ds.coords)

    # Save final output
    final_ds.to_netcdf(output_file)

    print("Final longitude ordering:", combined_ds['lon'].values)
    print(f"Combined NetCDF file saved to {output_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--workdir', type=str, required=True, help='Path to the working directory')
    parser.add_argument('--year', type=int, required=True, help='Year of the data to combine')
    parser.add_argument('--nametag', type=str, required=True, help='experiment name')
    args = parser.parse_args()

    file_pattern = f"{args.workdir}/{args.nametag}/{args.year}.era5-fluxes.nc.*"
    output_file = f"{args.workdir}/{args.year}.{args.nametag}-cleansky-spectra.nc"

    input_files = sorted(glob.glob(file_pattern))
    combine_netcdf_files(input_files, output_file)
