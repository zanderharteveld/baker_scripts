#######################################################
# this script looks at pdbs and checks if certain     #
# h-bond is present, if it is, the pdb is stored      #
# with others in folder, then best hbond score is     #
# copied to new folder                                #
# don't recommend using for without extensive editing #
# not good for general use                            #
#######################################################

import glob
import pyrosetta as py
import os
import sys
import shutil

py.init()
res1 = 22
res2 = 47
hb_dist = 4.0 
sf = py.get_fa_scorefxn()
final_folder = 'best_hb'
if not os.path.exists(final_folder):
        os.mkdir(final_folder)
score_manager = py.rosetta.core.scoring.ScoreTypeManager()
hb_score_term = score_manager.score_type_from_name("hbond_lr_bb")
for folder in glob.glob('*/'):
	max_hb_lr = 0
        max_hb_pdb = ''
        for pdb_file in glob.glob(folder + '*.pdb'):
#for pdb_file in glob.glob('11_0075/*.pdb'):
                pose = py.pose_from_pdb(pdb_file)
                N_res1 = pose.residue(res1).xyz("N")
                O_res1 = pose.residue(res1).xyz("O")
                N_res2 = pose.residue(res2).xyz("N")
                O_res2 = pose.residue(res2).xyz("O")
                first_hb_vector = N_res1 - O_res2
                second_hb_vector = N_res2 - O_res1
                if (first_hb_vector.norm() > hb_dist) and (second_hb_vector.norm() > hb_dist):
			print "Splited!"
		else:
			print "Not splited!"
                        if sf.score_by_scoretype(pose,hb_score_term,True) < max_hb_lr:
                                max_hb_lr = sf.score_by_scoretype(pose,hb_score_term,True)
                                max_hb_pdb = pdb_file
        split = max_hb_pdb.split('/')
        if len(split) > 1:
                file_name=split[0]+'_'+split[1]
                new_loc = "%s/%s" %(final_folder,file_name)
                shutil.copy(max_hb_pdb, new_loc)
