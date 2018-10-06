#!/bin/bash
#SBATCH -p medium 
#SBATCH -n 18
#SBATCH -N 1 
#SBATCH --mem=32g
#SBATCH -o log

python ./RMSD_allxall_parallel.py $(ls ~/loop_grafting/func_loops/cluster_loops/loops/*.pdb)
