###############################################################
# this script takes the loop results from the loop library    #
# and gets the sequences of the mutated loops and rename the  #
# pdb files to match the new convention                       #
# currently, the name splitting/appending of loops is not     #
# general for every naming convention                         #
# in additon, the naming of sequences in fasta file is        #
# based on the numbering of loops for our sequence, may not   #
# work for every beta barrel (lines 33-40)                    #
# inputs are the basename and the name of output file         #
# containing the sequences                                    #
###############################################################

import glob
import os
import shutil
import sys
import subprocess

pdbs = glob.glob('*pdb')
basename = sys.argv[1]+'_'
output = sys.argv[2]+'.txt'
i = 1
for pdb in pdbs:
    num_basename = basename+str(i).zfill(2)
    fasta_command = 'python2.7 /home/reggiano/Rosetta/main/source/scripts/python/public/pdb2fasta.py '+pdb
    fasta_out = subprocess.check_output(fasta_command.split(' '))
    fasta_seq = fasta_out.split()[-1]
    loops = list()
    edited_loops = list()
    name_split = pdb.strip('.pdb').split('_')
    loops.append(int(name_split[8]))
    if loops[-1] < 20:
       edited_loops.append(1)
    elif loops[-1] < 50:
       edited_loops.append(3)
    elif loops[-1] < 80:
       edited_loops.append(5)
    else:
       edited_loops.append(7)
    loops.append(int(name_split[9]))
    loops.append(int(name_split[12]))
    if loops[-1] < 20:
       edited_loops.append(1)
    elif loops[-1] < 50:
       edited_loops.append(3)
    elif loops[-1] < 80:
       edited_loops.append(5)
    else:
       edited_loops.append(7)
    loops.append(int(name_split[13]))
    #print(loops)
    os.rename(pdb,num_basename+'.pdb')
    with open(output, 'a') as f_out:
        f_out.write('>'+num_basename+'_loop_'+str(edited_loops[0])+'\n')
        f_out.write(fasta_seq[loops[0]-1:loops[1]]+'\n')
        f_out.write('>'+num_basename+'_loop_'+str(edited_loops[1])+'\n')
        f_out.write(fasta_seq[loops[2]-1:loops[3]]+'\n')
    f_out.close()
    i += 1
