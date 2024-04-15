#!/usr/bin/bash
#SBATCH --job-name=tdse
#SBATCH --nodes=2                     # number of nodes to use
#SBATCH --overcommit
#SBATCH --cpus-per-task=128           # number of CPUs PER NODE
#SBATCH --ntasks=4                   # number of total calculation processes
#SBATCH --mem=65536                   # typical simu takes 1.5GB during eigensolving phase
#SBATCH --output=logs/%j.log
#SBATCH --error=logs/%j.err
#SBATCH --export=ALL
#SBATCH --partition=cpu

pwd; hostname; date

echo "Running CPP TDSE on $SLURM_CPUS_ON_NODE CPU cores"

source /home/jmann/TDSESolveLinux/.venv/bin/activate

mpirun -x OMP_NUM_THREADS=64 -np 2 python3 -m mpi4py.futures $1