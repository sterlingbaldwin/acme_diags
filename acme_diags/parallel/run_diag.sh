#!/bin/bash
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -D $HOME/slurm_output/small/
#SBATCH -J plotset
#SBATCH -o set5_driver.o%j

source activate 2.8
ps5_commands.sh