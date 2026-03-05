module load cray-python/3.11.5

# Define experiment and year arrays
exps=(control PI co2_PI ch4_PI n2o_PI cfc12eq_PI hfc134aeq_PI co2_2xPI co2_4xPI)
years=($(seq 2001 2024))

# Calculate (exp, year) index from SLURM_ARRAY_TASK_ID
exp_index=$((SLURM_ARRAY_TASK_ID / ${#years[@]}))
year_index=$((SLURM_ARRAY_TASK_ID % ${#years[@]}))

workdir="/gpfs/f5/gfdl_m/scratch/Jing.Feng/work/grtcode-era5"
script="/ncrc/home1/Jing.Feng/scripts/GRTworkflow/combiner.py"


outdir="$workdir/$exp"

for exp in "${exps[@]}"; do
    for year in "${years[@]}"; do
        outdir="$workdir/$exp"
        echo "[$(date)] Running srun for $exp $year"
        
        srun python3 "$script" --workdir "$workdir" --year "$year" --nametag "$exp"
        
        if [ $? -eq 0 ]; then
            rm -f "$outdir/${year}"*.nc
            echo "[$(date)] Success: Removed $outdir/${year}*.nc"
        else
            echo "[$(date)] Failure: Skipped removal for $exp $year"
        fi
    done
done