#!/bin/bash -l
# created: Jun 22, 2018 3:37 PM
# author: minunnof
#SBATCH --job-name=calP9
#SBATCH --partition=small
#SBATCH --account=project_2000994
#SBATCH --output=sbaOut/output_%j.txt
#SBATCH --error=sbaOut/errors_%j.txt
#SBATCH --time=40:40:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10000
#SBATCH --mail-type=END
#SBATCH --mail-user=ritika.srinet@helsinki.fi

# set the number of threads based on --cpus-per-task
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
module load r-env-singularity
# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
    sed -i '/OMP_NUM_THREADS/d' ~/.Renviron
fi
# Specify a temp folder path
echo "TMPDIR=/scratch/project_2000994/Rtmp_heap" >> ~/.Renviron
# Run the R script
srun singularity_wrapper exec Rscript --no-save bcP9.R
# Correct permissions
chmod -Rf a+rwx /scratch/project_2000994/calibrations/srinet/newCal/chains || true

