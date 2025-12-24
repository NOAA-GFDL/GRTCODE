#!/bin/bash -e

# Get the directory where the environment script lives and the current time.
workflow_home="$( dirname -- "$( readlink -f -- "$0"; )"; )"
current_time="$( date +%Y-%m-%d-%H:%M:%S )"
grtcode_repository="$workflow_home/../grtcode"

# Variables needed to set reasonable defaults.
grtcode_input_data="$workflow_home/../grtcode/grtcode-data"

# Default arguments.
account="$( sacctmgr --noheader list user $USER format=DefaultAccount%6 )"
cfc11="$grtcode_input_data/cfc_cross_sections/CFC-11_absorption_cross_sections.csv"
cfc12="$grtcode_input_data/cfc_cross_sections/CFC-12_absorption_cross_sections.csv"
cluster="c5"
hitran_par_data="$grtcode_input_data/HITRAN_files/hitran2016.par"
h2o_continuum="$grtcode_input_data/water_vapor_continuum"
nodes="10"
n2_n2="$grtcode_input_data/collision_induced_absorption/N2-N2.csv"
o2_n2="$grtcode_input_data/collision_induced_absorption/O2-N2.csv"
o2_o2="$grtcode_input_data/collision_induced_absorption/O2-O2.csv"
o3_continuum="$grtcode_input_data/ozone_continuum/ozone_continuum.csv"
partition="batch"
queue="normal"
rfmip_irf_data="$grtcode_input_data/rfmip-irf/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc"
solar_flux_data="$grtcode_input_data/solar_flux/solar_flux.csv"
time_in_minutes="15"
threads="128"
work_directory="/gpfs/f5/gfdl_m/scratch/$USER/work/grtcode-work-directory-exp$experiment"
log="$work_directory/run-rfmip-irf-logfile-$current_time"

# Handle command line arguments.
argument_list="$0 [-h|--help] [-p <path to grtcode repository>] experiment"
counter=0
while [[ $# -gt 0 ]]; do
  argument="$1"
  case $argument in
    -h|--help)
      echo "$argument_list"
      echo "\nPositional arguments:"
      echo "experiment:    RFMIP-IRF experiment to run."
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
        experiment="$argument"
        if [ $experiment -lt 1 ] || [ $experiment -gt 18 ]; then
          echo "Error: experiment argument must be 1 <= experiment <= 18"
          exit 1
        fi
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
executable="$grtcode_repository/build/rfmip-irf"

# Create the runscript.
runscript="rfmip-irf-runscript.experiment-$experiment.$current_time"
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
    $rfmip_irf_data \\
    $(( experiment-1 )) \\
    -H2O -CO2 -O3 -N2O -CH4 -CO -O2 \\
    -h2o-ctm $h2o_continuum \\
    -o3-ctm $o3_continuum \\
    -CFC-11-eq $cfc11 \\
    -CFC-12 $cfc12 \\
    -N2-N2 $n2_n2 \\
    -O2-N2 $o2_n2 \\
    -O2-O2 $o2_o2 \\
    -integrated \\
    -flux-at-level 30 \\
    -x \$(( (\$index-1)*$nodes )) \\
    -X \$(( (\$index*$nodes)-1 )) \\
    -r-lw 0.1 \\
    -o $experiment.rfmip-irf-fluxes.forcing_index2.nc.\$index \\
    |& tee $experiment.rfmip-irf.log.\$index &

  if [ "\$?" -ne "0" ]; then
    echo "GRTCODE failed. See $work_directory/$experiment.rfmip-irf.log.\$index" > $log
    exit 1
  fi
done
wait
echo "Finished RFMIP-IRF experiment $experiment.  Combining the output." > $log

# Combine the output.
echo "Combining output files $work_directory/$experiment.rfmip-irf-fluxes.nc.xxxx" > $log
srun -v --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 \\
  $workflow_home/grtcode-results-combiner \\
  $work_directory/rfmip-irf.experiment-$experiment.nc \\
  $work_directory/$experiment.rfmip-irf-fluxes.forcing_index2.nc.1 \\
  10 \\
  |& tee log-combine
if [ "\$?" -ne "0" ]; then
  echo "Combining output files failed.  See ${work_directory}/log-combine" > $log
  exit 1
fi
echo "Output located at: $work_directory/rfmip-irf.experiment-$experiment.nc" > $log

# Delete temporary files.
for index in {1..$nodes}; do
  rm ${work_directory}/$experiment.rfmip-irf-fluxes.forcing_index2.nc.\$index
done
EOF

# Make the runscript executable.
chmod 755 $runscript

# Submit the runscript to slurm.
sbatch $runscript

# Print the results to the terminal too.
echo "See output from the job in $log once it finishes."
