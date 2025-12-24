#!/bin/bash -e

# Get the directory where the environment script lives and the current time.
workflow_home="$( dirname -- "$( readlink -f -- "$0"; )"; )"
current_time="$( date +%Y-%m-%d-%H:%M:%S )"
grtcode_repository="$workflow_home/../grtcode"

# Variables needed to set reasonable defaults.
grtcode_input_data="/ncrc/home1/Jing.Feng/scripts/grtcode/grtcode-data"

# Default arguments.
account="$( sacctmgr --noheader list user $USER format=DefaultAccount%6 )"
beta_distribution_data="/gpfs/f5/gfdl_m/scratch/Jing.Feng/line-by-line/run/input/clouds/beta_distribution.nc"
hfc134="$grtcode_input_data/cfc_cross_sections/HFC-134a_absorption_cross_sections.csv"
cfc12="$grtcode_input_data/cfc_cross_sections/CFC-12_absorption_cross_sections.csv"
cluster="c5"
era5_data="/gpfs/f5/gfdl_m/scratch/Jing.Feng/line-by-line/run/era5_test"
ghg_path="/gpfs/f5/gfdl_m/world-shared/Jing.Feng/GHG/annual_mean_gas_data.nc"
hitran_par_data="$grtcode_input_data/HITRAN_files/hitran2016.par"
h2o_continuum="$grtcode_input_data/water_vapor_continuum"
ice_cloud_parameterization_data="/gpfs/f5/gfdl_m/world-shared/Jing.Feng/cloud_optics/lbl_pade_ice_lw_solid_column_severlyroughen_gamma_aeq1_thick.nc"
liquid_cloud_parameterization_data="/gpfs/f5/gfdl_m/world-shared/Jing.Feng/cloud_optics/lbl_pade_liq_lw_mie_gamma_aeq12_thick.nc"
nodes="1"
n2_n2="$grtcode_input_data/collision_induced_absorption/N2-N2.csv"
o2_n2="$grtcode_input_data/collision_induced_absorption/O2-N2.csv"
o2_o2="$grtcode_input_data/collision_induced_absorption/O2-O2.csv"
o3_continuum="$grtcode_input_data/ozone_continuum/ozone_continuum.csv"
partition="batch"
queue="normal"
solar_flux_data="$grtcode_input_data/solar_flux/solar_flux.csv"
time_in_minutes="20"
threads="128"

# Handle command line arguments.
argument_list="$0 [-h|--help] [-p <path to grtcode repository>] year"
counter=0
while [[ $# -gt 0 ]]; do
  argument="$1"
  case $argument in
    -h|--help)
      echo "$argument_list"
      echo "\nPositional arguments:"
      echo "year:          ERA5 year to run."
      echo "\nOptional arguments:"
      echo "-h, --help:    Print this message."
      echo "-p <path>:     Path to a grtcode repository."
      exit 0
    ;;
    -p)
      shift
      if [ "$#" -gt 0 ]; then
        grtcode_repository="$(cd "$1" && pwd)"
        shift
      else
        echo "Error: please specify the path to the grtcode repository."
        echo "Usage: $argument_list\n"
        exit 1
      fi
    ;;
    *)
      counter=$((counter+1))
      if [ $counter -eq 1 ]; then
        year="$argument"
      else
        echo "Error: too many arguments."
        echo "Usage: $argument_list"
        exit 1
      fi
      shift
    ;;
  esac
done
if [ $counter -ne 1 ]; then
  echo "Error: missing required argument."
  echo "Usage: $argument_list"
  exit 1
fi
executable="$grtcode_repository/build/era5"
work_directory="/gpfs/f5/gfdl_m/scratch/$USER/work/grtcode-era5-$year"
log="$work_directory/run-era5-logfile-$current_time"

# Create the runscript.
runscript="era5-runscript.$year.$current_time"
cat <<EOF >> $runscript
#!/bin/bash -ex

#SBATCH --account=$account
#SBATCH --partition=$partition
#SBATCH --clusters=$cluster
#SBATCH --qos=$queue
#SBATCH --nodes=$nodes
#SBATCH --time=$time_in_minutes

# Load the software environment.
source $workflow_home/environment.sh

# Create a work directory and go there.
mkdir -p $work_directory
cd $work_directory
touch $log

# Run GRTCODE.
export OMP_NUM_THREADS=$threads
for index in {1..$nodes}; do
  srun -v --nodes=1 --ntasks-per-node=1 --cpus-per-task=128 \\
    $executable \\
    $hitran_par_data \\
    $solar_flux_data \\
    $era5_data/$year-era5.nc \\
    $greenhouse_gas_data \\
    -clean \\
    -year $year \\
    -H2O -CO2 -O3 -N2O -CH4 -CO -O2 \\
    -h2o-ctm $h2o_continuum \\
    -o3-ctm $o3_continuum \\
    -HFC-134a-eq $hfc134 \\
    -CFC-12-eq $cfc12 \\
    -N2-N2 $n2_n2 \\
    -O2-N2 $o2_n2 \\
    -O2-O2 $o2_o2 \\
    -x \$(( (\$index-1)*(144/$nodes) )) \\
    -X \$(( (\$index*(144/$nodes)-1) )) \\
    -o $year.era5-fluxes.nc.\$index \\
    -beta-path $beta_distribution_data \\
    -ice-path $ice_cloud_parameterization_data \\
    -liquid-path $liquid_cloud_parameterization_data \\
    -flux-at-level 12 \\
    -r-lw 1\\
    -r-sw 10\\
    -ghg_start_year 1750 \\
    |& tee $year.era5.log.\$index &

  if [ "\$?" -ne "0" ]; then
    echo "GRTCODE failed. See $work_directory/$year.era5.log.\$index" > $log
    exit 1
  fi
done
wait
echo "Finished ERA5 $year.  Combining the output." > $log

# Combine the output.
echo "Combining output files $work_directory/$year.era5-fluxes.nc.xxxx" > $log

# Run the Python script
srun -v --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 \
  python3 //ncrc/home1/Jing.Feng/scripts/GRTworkflow/combiner.py --workdir "$work_directory" --year "$year" \
  |& tee "$log"

if [ "\$?" -ne "0" ]; then
  echo "Combining output files failed.  See ${work_directory}/log-combine" > $log
  exit 1
fi
echo "Output located at: $work_directory/era5-$year.nc" > $log

# Delete temporary files.
#for index in {1..$nodes}; do
#  rm ${work_directory}/$year.era5-fluxes.nc.\$index
#done
EOF

# Make the runscript executable.
chmod 755 $runscript

# Submit the runscript to slurm.
sbatch $runscript

# Print the results to the terminal too.
echo "See output from the job in $log once it finishes."
