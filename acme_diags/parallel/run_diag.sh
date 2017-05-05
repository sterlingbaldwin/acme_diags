#!/bin/bash
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -D /opt/nfs/mcenerney1/slurm_output/small/
#SBATCH -J plotset
#SBATCH -o set5_driver.o%j

module use /usr/common/contrib/acme/modulefiles
module load uvcdat/batch

ps5_commands.sh