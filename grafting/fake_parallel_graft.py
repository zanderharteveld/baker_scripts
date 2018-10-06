#######################################################
# this script will submit all possible combinations   #
# of loops based on the loop positions given and the  #
# loops contained in the nested directory             #
# jobs are submitted as array task, run.sh and        #
# loop_grafter_renum_min.py needs                     #
# to be in folder in order for this to run            #
# really recommend clustering your loops by RMSD      #
# before running this script                          #
# example command call is shown in:                   #
# /home/reggiano/scripts/grafting/submit_all_cmd      #
#######################################################


import os
import shutil
import glob
import sys
import argparse
import fileinput
from random import shuffle
from math import ceil


def main():
    args = parseargs()
    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)
    all_commands = list()
    all_loop_dirs = glob.glob(args.nested_dir+'/*')
    for i in range(0,len(args.residue_ranges)-1):
        for j in range(i+1,len(args.residue_ranges)):
            loop1 = args.residue_ranges[i].strip("[]").split(",")[0]+' '+args.residue_ranges[i].strip("[]").split(",")[1]
            loop2 = args.residue_ranges[j].strip("[]").split(",")[0]+' '+args.residue_ranges[j].strip("[]").split(",")[1]
            res1_file_str = loop1.split(' ')[0]+'_'+loop1.split(' ')[1]
            res2_file_str = loop2.split(' ')[0]+'_'+loop2.split(' ')[1]
            #all combinations of loop directories
            for m in range(0,len(all_loop_dirs)-1):
                for n in range(m+1,len(all_loop_dirs)):
                    if m != n:
                        ldir1 = all_loop_dirs[m]
                        ldir2 = all_loop_dirs[n]

                        #first submit is with loopdir1,loopdir2 in loop1, loop2 positions
                        outer_dir = args.out_dir+'/'+ldir1.split('/')[-1]+'_'+res1_file_str
                        if not os.path.exists(outer_dir):
                            os.mkdir(outer_dir)
                        output_dir = outer_dir+'/'+ldir1.split('/')[-1]+'_'+res1_file_str+'_'+ldir2.split('/')[-1]+'_'+res2_file_str
                        command_lines = store_run_commands(args, ldir1, ldir2, loop1, loop2, output_dir)
                        all_commands.append(command_lines)

                        #second submit is with loopdir2, loopdir1 in loop1, loop2 positions
                        outer_dir = args.out_dir+'/'+ldir2.split('/')[-1]+'_'+res1_file_str
                        if not os.path.exists(outer_dir):
                            os.mkdir(outer_dir)
                        output_dir = outer_dir+'/'+ldir2.split('/')[-1]+'_'+res1_file_str+'_'+ldir1.split('/')[-1]+'_'+res2_file_str
                        command_lines = store_run_commands(args, ldir2, ldir1, loop1, loop2, output_dir)
                        all_commands.append(command_lines)
    
    #shuffle commands to randomize order
    shuffle(all_commands)

    #prepare array tasks, 10000 in each list
    num_tasks_lists = int(ceil(len(all_commands)/10000))+1
    print('total # of commands: %d' %len(all_commands))
    print('total # of task lists: %d' %num_tasks_lists)
    j = 0
    file_name = "graft_tasks_"
    for i in range(1,num_tasks_lists+1):
        task_file = file_name+str(i)
        with open(task_file, 'w') as f:
            first_one = True
            while ( ( j % 10000 != 0) or first_one ) and ( j < len(all_commands) ):
                if first_one:
                    first_one = False
                f.write(all_commands[j])
                j += 1
        f.close()
        first_one = True

    #submit each array batch job
    task_files = glob.glob(file_name+'*')
    for task_list in task_files:     
        update_run(args, task_list)
        submit_cmd = 'sbatch -a 1-$(cat '+task_list+'| wc -l) run.sh'
        #os.system('cat run.sh')
        print(submit_cmd)
        #os.system(submit_cmd)
                      
def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', type=str, help='input pdb file to graft onto')
    parser.add_argument('-nd', '--nested_dir', type=str, help='directory containing all loops you want to use')
    parser.add_argument('-r', '--residue_ranges', type=str, nargs='+', help='all the loop anchors from orig scaffold you want to graft on.'
                                                       'list like so, C-terminal loops first: [90,100] [24,30] ... [1,5]')
    parser.add_argument('-sc', '--shape_comp', type=str, help='are you optimizing for ligand shape complementarity?')
    parser.add_argument('-l', '--ligand', type=str, help='if neccessary, the ligand params file')
    parser.add_argument('-p', '--queue', type=str, default="short", help='queue you want to submit to')
    parser.add_argument('-o', '--out_dir', type=str, default='grafted_loops',help='name of directory containing all output_dir')
    args = parser.parse_args()
    return args

#make a single list of all important lines for single call
def store_run_commands(args, loop_dir1, loop_dir2, anchor1, anchor2, out_dir):
    python_call = "python loop_grafter_renum_min.py -i %s -d %s %s -o %s -a1 %s -a2 %s -sc %s -l %s \n" \
                            %(args.input_file, loop_dir1, loop_dir2, out_dir, anchor1, anchor2, args.shape_comp, args.ligand) 
    return python_call

#updates run script with new arguments
def update_run(args, task_list):
    for line in fileinput.input('run.sh', inplace=1):
        if '#SBATCH -p' in line:
            line = '#SBATCH -p '+args.queue +'\n'
        if '#SBATCH -o' in line:
            line = '#SBATCH -o graft_list_'+task_list.split('_')[-1]+'.log \n'
        if 'CMD=$' in line:
            line = 'CMD=$(head -n $SLURM_ARRAY_TASK_ID '+task_list+' | tail -1) \n'
        sys.stdout.write(line)

main()
