######################################################
# this script looks up the abego type of turn from   #
# residues marked by segment lookup                  #
# then it gets important residues for design based   #
# on frequency seen in pdb                           #
# that step needs to be optimized, very hacky right  #
# now                                                #
# other scripts in this folder must be in the        #
# cwd when this is called                            #
###################################################### 

import pyrosetta as py
import glob
import os
import sys
import subprocess as sp
import numpy as np
import shutil
py.init()

def main():
    pdbs = glob.glob('*pdb')
    abego_man = py.rosetta.core.sequence.ABEGOManager()
    dssp = py.rosetta.protocols.moves.DsspMover()
    
    #get aa dictionary for later
    aa_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    aa_dict={}

    if not os.path.exists('abego_design'):
        os.mkdir('abego_design')

    for i,aa in enumerate(aa_list):
        aa_dict[aa] = i
        
    for pdb in pdbs:
        seg_look_residues = get_residues(pdb)
        pose = py.pose_from_pdb(pdb)
        dssp.apply(pose)
        ss = str(pose.secstruct())
        loop_abego = ''
        loop_res = list()
        full_abego = py.rosetta.core.sequence.ABEGOManager.get_symbols(abego_man, pose)
        #print(full_abego)
        #keep track of residues involved in abego and get abego determination
        copy = False
	for res in seg_look_residues:
            if ((ss[res-2] == 'L' or ss[res-2] == 'H') and ss[res-1] == 'E' and ss[res] == 'E' and copy):
                loop_res.append(res)
#                print('here 3 '+str(res))
#                print(copy)
	        break
            if copy:
                loop_abego+=str(full_abego[res])
                loop_res.append(res)
#                print('here 1 '+str(res))
#                print(copy)
            if ((ss[res] == 'L' or ss[res] == 'H') and ss[res-1] == 'E' and ss[res-2] == 'E' and not copy):
                copy = True
	        loop_res.append(res-1)
                loop_res.append(res)
                loop_abego+=str(full_abego[res])
#                print(copy)
#         	print('here 2 '+str(res))

        #call and keep output from abego python script
        command = 'python read_stats_HH_cap-nograph.py -data abego_seq_stats_PDB_HH-EH-HE-EE-wcap.out -connection EE -abego %s' %loop_abego
        p = sp.Popen(command, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, close_fds=True)
        #stdout = p.communicate()[0]
        abego_list = list()
        abego_freq = 0
        #get output into list
        line_num = 0
        #print(loop_abego)
        #print(loop_res)
        while True:
            line = p.stdout.readline()
            #print(line)
            if line != '':
                line_num += 1
            	if line_num >= 4:
                	line = line.split()
                	abego_list.append(line[1])
                	abego_freq += 1.0
	    else: 
                break
                
        #calculate frequencies of amino acids at certain residues, get allowed_aa at that position
        npos = len(loop_abego)+2
        posfreq = np.zeros((npos,20))
        for abego in abego_list:
            #print(abego)
            if len(abego) > 0:
                for pos in range(npos):
                    aa = abego[pos]
                    aa_num = aa_dict.get(aa)
                    posfreq[pos][aa_num] += 1
        posfreq = posfreq/abego_freq
        #print(posfreq)
        allowed_aa = list()
        if len(abego_list) > 1: #checks that we are not restricting loop to a single structure
            for pos in range(npos):
                #mean = np.mean(posfreq,pos,)
                aa_at_pos = ''
                for aa_num in range(20):
                    aa = aa_list[aa_num]
                    if posfreq[pos][aa_num] >= 0.3:
                        aa_at_pos += aa
                allowed_aa.append(aa_at_pos)

        #copy pdbs to new folder, and write resfile in that new folder
        res_file = 'abego_design/'+pdb[:-4]+'.resfile'
        total_aa = 0
        for aa in allowed_aa:
                if aa != '':
                    total_aa += 1

        #checks that there was decent number output from abego lookup and that all aa spots aren't blank
        if len(abego_list) > 4 and total_aa > 0:
            new_loc = 'abego_design/'+pdb
            shutil.copy(pdb,new_loc)
            with open(res_file, 'w') as f_out:
               f_out.write('start\n')
               line_written = False
               for i,res in enumerate(loop_res):
                    if len(allowed_aa[i]) > 0:
                        if len(allowed_aa[i]) == 1:
                            aa_num = aa_dict.get(allowed_aa[i])
                        else:
                            max_freq = 0
                            for char in allowed_aa[i]:
                                aa_tmp = aa_dict.get(char)
                                if posfreq[i][aa_tmp] > max_freq:
                                    max_freq = posfreq[i][aa_num]
                                    aa_num = aa_tmp

                        #both 92 and 93 are conserved from the original design
                        #94 can be added in addition to conserved, if and only if the residue is important for the beta turn
                        if res == 92: 
                            line = str(res)+' A PIKAA L \n'
                            line_written = True
                            f_out.write(line)
                        if res == 93:
                            if not line_written:
                                line = '92 A PIKAA L \n'
                                f_out.write(line)
                                line_written = True
                            line = str(res)+' A PIKAA F \n'
                            f_out.write(line)
                        if res == 94:
                            if posfreq[i][aa_num] >= 0.8:
                                extra_aa = allowed_aa[i]
                            else:
                                extra_aa = ''
                            if not line_written:
                                line = '92 A PIKAA L \n93 A PIKAA F \n'
                                f_out.write(line)
                                line_written = True
                            line = str(res)+' A PIKAA K'+extra_aa+'\n'
                            f_out.write(line)
                        if res > 94 and not line_written:
                            line = '92 A PIKAA L \n93 A PIKAA F \n94 A PIKAA K\n'
                            f_out.write(line)
                            line = str(res)+' A PIKAA '+allowed_aa[i]+'\n'
                            f_out.write(line)
                            line_written = True
                        
                        #prevents letting the 94 from the new abego design through when it should not be
                        elif res > 94 and line_written: 
                            line = str(res)+' A PIKAA '+allowed_aa[i]+'\n'
                            f_out.write(line)
            f_out.close()

def get_residues(pdb_file):
    res=list()
    with open(pdb_file,'r') as pdb:
        for line in pdb:
            if 'REMARK' in line:
                line=line.split()
                res.append(int(line[2]))
    return(res)


main()
