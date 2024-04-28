#!/usr/bin/bash
#SBATCH --job-name=tdse
#SBATCH --nodes=1                     # number of nodes to use
#SBATCH --cpus-per-task=64           
#SBATCH --ntasks=1                   # number of total calculation processes
#SBATCH --mem=65536                   # typical simu takes 1.5GB during eigensolving phase
#SBATCH --output=logs/%j.log
#SBATCH --error=logs/%j.err
#SBATCH --export=ALL
#SBATCH --partition=cpu

pwd; hostname; date

echo "Running CPP TDSE on $SLURM_CPUS_ON_NODE CPU cores"

source /home/jmann/TDSESolveLinux/.venv/bin/activate

export OMP_NUM_THREADS=64

python3 -u $1