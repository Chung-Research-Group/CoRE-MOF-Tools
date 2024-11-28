#!/bin/bash

#SBATCH --partition=amd
<<<<<<< HEAD
#SBATCH -w SG13
=======
>>>>>>> 04b52f8f16bcd5c63e8428b2616e7fe32198d02e
#SBATCH --nodes=1
##SBATCH --ntasks=1
##SBATCH --ntasks-per-node=1
##SBATCH --ntasks-per-socket=1
##SBATCH --cpus-per-task=1
##SBATCH --time=10:00:00
#SBATCH --job-name=my_job

cd $SLURM_SUBMIT_DIR

export LD_LIBRARY_PATH=$RASPA_DIR/lib:$LD_LIBRARY_PATH
$RASPA_DIR/bin/simulate simulation.input
