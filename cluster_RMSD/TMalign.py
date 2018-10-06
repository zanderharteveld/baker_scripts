from subprocess import check_output
from os import path
from sys import argv

def get_RMSD(tmout):
    lines = tmout.splitlines()
    rmsd_line = lines[16]
    rmsd = rmsd_line.split()[4].split(',')[0]
    return float(rmsd)

def get_tmscore(tmout):
    lines = tmout.splitlines()
    tmscore1_line = lines[17]
    tmscore2_line = lines[18]
    tmscore1 = float(tmscore1_line.split()[1])
    tmscore2 = float(tmscore2_line.split()[1])
    tmscore = ( tmscore1 + tmscore2 )/2
    return tmscore

def TMalign(str1,str2,**kwargs):
    to_exec = "/work/basantab/bin/TMalign"
    if 'exe' in kwargs: # get the TMalign executable path from optional args
        to_exec = kwargs['exe']
    TM_out = check_output([to_exec, str1, str2]).decode('ascii')
    return TM_out
        
if __name__ == '__main__':
    # test with:
    # python TMalign.py ~/NTF2_project/20151125_3x4_len1_H3x2_seq/only_freq_ABEGOs/hyak_run_20151208/20151208_sequence_on_short_big_pocket/j0/0/1573_L1H18L2H7L5E4L2E4L1H12L1H8L1E11L2E10L2E12L2E10L1_41_4_input_0006_1filtered_0001.pdb ~/NTF2_project/20151125_3x4_len1_H3x2_seq/only_freq_ABEGOs/hyak_run_20151208/20151208_sequence_on_short_big_pocket/j0/0/1573_L1H18L2H7L5E4L2E4L1H12L1H8L1E11L2E10L2E12L2E10L1_41_4_input_0006_1filtered_0001_0015.pdb ~basantab/bin/TMalign
    str1 = argv[1]
    str2 = argv[2]
    path = argv[3]
    out = TMalign(str1,str2)
    out = TMalign(str1,str2,exe=path)
    print(get_RMSD(out))
