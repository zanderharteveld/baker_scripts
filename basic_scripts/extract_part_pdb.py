####################################
# this script extracts part of pdb #
# and saves as a new pdb           #
# variables need to be edited      #
####################################

import os
import glob

loop_name = "L3A"
loop_start = 44
loop_end = 48
for pdb in glob.glob("*.pdb"):
	loop_number = int(os.path.splitext(pdb)[0])
	pdb_out = loop_name + "_" + "%04d.pdb" %(loop_number)
	with open(pdb, 'r') as f_in:
		with open(pdb_out, 'w') as f_out:
			for line in f_in:
					if line[0:4] == "ATOM":
						if (int(line.split()[5]) >= loop_start) and (int(line.split()[5]) <= loop_end) :
							f_out.write(line)
			
			
