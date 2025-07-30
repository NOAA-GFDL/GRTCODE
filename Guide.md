# GRTcode

GRTcode is a radiative transfer modeling tool used for IRF (Instantaneous Radiative Forcing) calculations and climate diagnostics. This repository provides the core codebase, workflows for running experiments, and example datasets.

## Clone the Repository

```bash
git clone https://github.com/fengzydy/GRTCODE.git
cd GRTCODE
```
## Build Instructions (on GAEA)

To compile the GRTcode executables:

```bash
./GRTworkflow/build.sh -p ./
```

> 💡 If you omit the optional `-p <path>` argument, the script will only report the status of the pre-built executables.


## Prepare Input Spectral Data

```bash
./download-test-data
unzip grtcode-data
```

This will download and extract example spectral datasets needed for test runs.

## Run Examples

### Run RFMIP-IRF 100-Column Experiments

```bash
./GRTworkflow/run-rfmip-irf.sh -p <path_to_GRTcode> <experiment_number>
```

Example (to run experiment 1):

```bash
./GRTworkflow/run-rfmip-irf.sh -p ./ 1
```

### Run One Year of ERA5 Reanalysis

```bash
./GRTworkflow/run-era5.sh -p ./ <year>
```

Example:

```bash
./GRTworkflow/run-era5.sh -p ./ 2010
```

---

### Run Multiple ERA5 Experiments

To run a series of experiments over multiple years:

1. Open and edit the submission script:

```bash
vim ./GRTworkflow/submit.sh
```

2. Modify the paths for `runscript` and executables accordingly.

3. Example loop inside `submit.sh`:

```bash
for year in $(seq 2010 2010); do
    for exp in "${exps[@]}"; do
        cmd="$runscript -p /ncrc/home1/Jing.Feng/scripts/grtcode $year $era5_data $data_dir/${exp}.nc ${exp}_noh2o"
        echo "Running: $cmd"
        eval "$cmd"
    done
done
```

To submit the batch job:

```bash
./GRTworkflow/submit.sh
```

## Download ERA5 Input Data

Navigate to the `era5-tools` directory and run:

```bash
cd era5-tools
python3 download_2021to2025.py
```

> Edit the script to modify years and output directory as needed.


## 🧰 Python Environment Setup 

```bash
module load cray-python/3.11.5
pip install xarray glob os argparse netCDF4
```

