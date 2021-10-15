#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p serial
#SBATCH -o result.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtg374@nyu.edu
#SBATCH -t 3-00:00

module purge
module load matlab/2018a
Cmd=$(printf '%s' "$1" )
matlab -nodisplay -r "$Cmd;quit"
