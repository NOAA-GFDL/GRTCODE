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

        # Filter variables starting with 'rl' and convert float64 -> float32
        rl_vars = {}
        for var in ds.data_vars:
            if var.startswith('rl'):
                data = ds[var]
                # Replace NaN _FillValue with -9999.0
                if '_FillValue' in data.attrs and np.isnan(data.attrs['_FillValue']):
                    data.attrs['_FillValue'] = -9999.0
                # Convert dtype if needed
                if data.dtype == 'float64':
                    data = data.astype('float32')
                rl_vars[var] = data

        # Also cast coords if float64, and sanitize _FillValue
        coords = {}
        for name, coord in ds.coords.items():
            if coord.dtype == 'float64':
                coord = coord.astype('float32')
            if '_FillValue' in coord.attrs and np.isnan(coord.attrs['_FillValue']):
                coord.attrs['_FillValue'] = -9999.0
            coords[name] = coord

        rl_ds = xr.Dataset(rl_vars, coords=coords)
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

    # Reorder dimensions if all present
    desired_order = ('time', 'lw_wavenumber', 'lat', 'lon')
    current_dims = list(combined_ds.dims)
    
    try:
        if set(desired_order).issubset(set(current_dims)):
            # Add any remaining dims to ensure all are present
            remaining_dims = [dim for dim in current_dims if dim not in desired_order]
            full_order = list(desired_order) + remaining_dims
            # Transpose only if full_order is a permutation of current_dims
            if set(full_order) == set(current_dims):
                combined_ds = combined_ds.transpose(*full_order)
            else:
                print(f"Skipping transpose due to mismatch: {full_order} vs {current_dims}")
        else:
            print(f"Skipping transpose: dataset dims = {current_dims}, expected at least = {desired_order}")
    except Exception as e:
        print(f"Transpose failed with error: {e}")

    combined_ds.to_netcdf(output_file)
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
