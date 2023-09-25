#!/bin/bash
#SBATCH --job-name=15_E_DE_batch
#SBATCH --time=24:00:00
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G
#SBATCH --mail-type=end
#SBATCH --mail-user=jmitch81@jhmi.edu
#SBATCH -o ./report/slurm-%A_%a.out

module load seurat/4.1.1
module list

R_SCRIPT="scripts/15_E_Differential_Expression_Batch_Control.R"

R CMD BATCH $R_SCRIPT

echo "Finished with job $SLURM_JOBID"
