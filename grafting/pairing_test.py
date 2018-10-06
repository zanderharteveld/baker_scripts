##################################################
# from anatassia, general script for checking if #
# certain hbond is present in pdbs               #
##################################################

import glob
import pyrosetta as py
import os
import sys
import shutil

py.init()


for folder in glob.glob('*/'):
	for pdb_file in glob.glob(folder + '*.pdb'):
#for pdb_file in glob.glob('11_0075/*.pdb'):
                pose = py.pose_from_pdb(pdb_file)
                hbond_set = py.rosetta.core.scoring.hbonds.HBondSet()
                pose.update_residue_neighbors()
                py.rosetta.core.scoring.hbonds.fill_hbond_set(pose, False, hbond_set)
                O92_N18 = False
                O18_N92 = False
#		print hbond_set.nhbonds()
#		print hbond_set.show(pose)

                for i in range(1, hbond_set.nhbonds()+1):
                        if hbond_set.hbond(i).acc_res() == 22 and hbond_set.hbond(i).acc_atm_is_backbone() == True and hbond_set.hbond(i).don_res() == 47 and hbond_set.hbond(i).don_hatm_is_backbone() == True:
				O18_N92 = True
			elif hbond_set.hbond(i).don_res() == 22 and hbond_set.hbond(i).don_hatm_is_backbone() == True and hbond_set.hbond(i).acc_res() == 47 and hbond_set.hbond(i).acc_atm_is_backbone() == True:
                                O92_N18 = True

		if O92_N18 == False and O18_N92 == False:
			print "Splited!"
		else:
			print "Not splited!"
                        split = pdb_file.split('/')
                        file_name=split[0]+'_'+split[1]
                        new_loc = "not_splited/%s" %(file_name)
                        shutil.copy(pdb_file, new_loc)
