########################################
# this script extracts the pssm        #
# from the end of the pdb after        #
# segment lookup                       #
########################################

import os
import glob

for pdb in glob.glob("*.pdb"):
	begin = False
	end = False
	pssm = '%s.pssm' %(pdb[:-4])
	with open(pdb, 'r') as f_in:
		with open(pssm, 'w') as f_out:
			for line in f_in:
				if "#BEGIN_SEGMENT_SEQUENCE_PROFILE_PSSM" in line:
					begin = True
				if begin:
					f_out.write(line)
				if "#END_SEGMENT_SEQUENCE_PROFILE_PSSM" in line:
					end = True
				if end:
					break
