#!/bin/bash

#SBATCH --partition=your_partition
#SBATCH -w your_node
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
##SBATCH --ntasks=1
##SBATCH --ntasks-per-node=1
##SBATCH --ntasks-per-socket=1
##SBATCH --cpus-per-task=1
##SBATCH --time=10:00:00
#SBATCH --job-name=my_job

cd $SLURM_SUBMIT_DIR

export LD_LIBRARY_PATH=$RASPA_DIR/lib:$LD_LIBRARY_PATH
$RASPA_DIR/bin/simulate simulation.input
