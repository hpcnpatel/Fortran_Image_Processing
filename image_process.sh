#!/bin/bash
# Job Name and Files (also --job-name)
#SBATCH -J image_process
#Output and error (also --output, --error):
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#Initial working directory (also --chdir):
#SBATCH -D ./
#Notification and type
#SBATCH --mail-type=END
#SBATCH --mail-user=patel@lrz.de
# Wall clock limit:
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#Setup of execution environment
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pr28fa
#constraints are optional
#--constraint="scratch&work"
#SBATCH --partition=micro
#Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=48
#Important
module load slurm_setup

#Run the program:
#mpiexec -n $SLURM_NTASKS echo "Hello world"
#mpiexec echo "Hello world"
#srun -n $SLURM_NTASKS echo "Hello world"

./inter parameter.dat
