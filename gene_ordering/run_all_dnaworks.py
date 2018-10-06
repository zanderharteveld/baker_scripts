############################################
# this script will run dnaworks for every  #
# pdb in folder, and also rename pdbs with #
# desired basename+numbers                 #
# inputs are the basname and the general   #
# dna template to copy from                #
############################################

import glob
import os
import shutil
import sys
import subprocess

pdbs = glob.glob('*pdb')
basename = sys.argv[1]+'_'
dna_input_template = sys.argv[2]
i = 1
for pdb in pdbs:
    num_basename = basename+str(i).zfill(4)
    fasta_command = 'python2.7 /home/reggiano/Rosetta/main/source/scripts/python/public/pdb2fasta.py '+pdb
    fasta_out = subprocess.check_output(fasta_command.split(' '))
    fasta_seq = fasta_out.split()[-1]
    os.rename(pdb,num_basename+'.pdb')
    with open(dna_input_template, 'r') as template:
        with open(num_basename+'.inp', 'w') as inp_out:
            for j, line in enumerate(template):
                if 'TITLE' in line:
                    line = 'TITLE '+num_basename+' \n'
                if 'LOGFILE' in line:
                    line = 'LOGFILE '+num_basename+'.log \n'
                if j == 25:
                    line = fasta_seq+'\n'
                inp_out.write(line)
            inp_out.close()
        template.close()
    dnaworks_cmd = '/home/strauch/local/DNAWorks/dnaworks '+num_basename+'.inp'
    #os.system(dnaworks_cmd)
    subprocess.call(dnaworks_cmd.split(' '))
    i += 1 
