#############################################
# takes all logs from dnaworks and extracts #
# the dnaseq and adds the seq for vector    #
# pcDB347                                   #
#############################################


import glob
import os

files = glob.glob('*log')
output ='all_dna_seq'
for file in files:
    seq_lines = list()
    with open(file, 'r') as f_in:
        for i,line in enumerate(f_in):
            if i >= 114 and i < 120:
                seq_lines.append(line)
    with open(output, 'a') as f_out:
        f_out.write('>'+file[:-4]+'\n')
        f_out.write('\n')
        f_out.write('GCGAAAACCTGTATTTTCAG')
        for line in seq_lines:
            seq = line.strip().split(' ')[-1]
            f_out.write(seq)
        f_out.write('TAATGACGGCTGCTAACAAAGCCCGA')
        f_out.write('\n')
        f_out.write('\n')
