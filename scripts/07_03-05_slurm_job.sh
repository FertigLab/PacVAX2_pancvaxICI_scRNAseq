#!/bin/bash
#SBATCH --job-name=03-05_cluster_DE
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=480G
#SBATCH --mail-type=end
#SBATCH --mail-user=jmitch81@jhmi.edu
#SBATCH -o ./report/slurm-%A_%a.out

module load seurat/4.1.1
module list

R_SCRIPT="scripts/07_03-05_Cluster_Differential_Expression.R"

R CMD BATCH $R_SCRIPT

echo "Finished with job $SLURM_JOBID"