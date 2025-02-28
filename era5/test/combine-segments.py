from argparse import ArgumentParser
from glob import glob
from shutil import copyfile

from netCDF4 import Dataset


def combine_segments(basename, output):
    files = glob(".".join([basename, "segment*"]))
    with Dataset(output, "w") as dataset:
        for i, path in enumerate(files):
            with Dataset(path, "r") as old:
                if i == 0:
                    for name, dim in old.dimensions.items():
                        if name == "lon":
                            size = old.getncattr("lon_global_size")
                        else:
                            size = dim.size
                        dataset.createDimension(name, size)
                    for name, v in old.variables.items():
                        dataset.createVariable(name, v.dtype, v.dimensions)
                        for att in v.ncattrs():
                            if att not in ["_FillValue",]:
                                dataset[name].setncattr(att, v.getncattr(att))
                for name, v in dataset.variables.items():
                    if "lon" in v.dimensions:
                        x = old.getncattr("lon_start")
                        X = old.getncattr("lon_stop") + 1
                        if name == "lon":
                            v[x:X] = old[name][...]
                        elif "wavenumber" in v.dimensions:
                            v[..., x:X, :] = old[name][...]
                        else:
                            v[..., x:X] = old[name][...]
                    else:
                        v[...] = old[name][...]


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("basename", help="base name for the segments.")
    parser.add_argument("output", help="name of output dataset.")
    args = parser.parse_args()
    combine_segments(args.basename, args.output)
