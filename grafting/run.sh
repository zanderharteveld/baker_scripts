#!/bin/bash
#SBATCH -p medium 
#SBATCH --mem=4g
#SBATCH -o graft_list_6.log 

CMD=$(head -n $SLURM_ARRAY_TASK_ID graft_tasks_6 | tail -1) 
source /software/miniconda3/envs/pyrosetta2/bin/activate pyrosetta3 
exec ${CMD}

