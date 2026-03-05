#!/bin/bash
#SBATCH --job-name=grt-combine
#SBATCH --account=gfdl_m
#SBATCH --partition=batch
#SBATCH --clusters=c5
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --time=00:05:00
#SBATCH --array=0-247
#SBATCH --output=logs/slurm-%A_%a.out
#SBATCH --error=logs/slurm-%A_%a.err

module load cray-python/3.11.5

# Define experiment and year arrays
exps=(control PI co2_PI ch4_PI n2o_PI cfc12eq_PI hfc134aeq_PI co2_2xPI co2_4xPI)
years=($(seq 2001 2024))

# Calculate (exp, year) index from SLURM_ARRAY_TASK_ID
exp_index=$((SLURM_ARRAY_TASK_ID / ${#years[@]}))
year_index=$((SLURM_ARRAY_TASK_ID % ${#years[@]}))

exp=${exps[$exp_index]}
year=${years[$year_index]}

workdir="/gpfs/f5/gfdl_m/scratch/Jing.Feng/work/grtcode-era5"
script="/ncrc/home1/Jing.Feng/scripts/GRTworkflow/combiner.py"
outdir="$workdir/$exp"

echo "[$(date)] Task $SLURM_ARRAY_TASK_ID — Running $exp $year"

# Run the combiner using srun
python3 "$script" --workdir "$workdir" --year "$year" --nametag "$exp"

# Only delete files if the command above succeeded
if [ $? -eq 0 ]; then
    rm -f "$outdir/${year}"*.nc
    echo "[$(date)] Success: Removed $outdir/${year}*.nc"
else
    echo "[$(date)] Failure: Skipped removal for $exp $year"
fi

