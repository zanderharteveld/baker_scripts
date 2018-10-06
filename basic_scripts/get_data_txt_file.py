##############################################
# this script is for CD spectra analysis     #
# extracts data from dat file into txt file  #
# readable by python                         #
##############################################

import glob
import os

for file in glob.glob("*dat"):
	new_file = file[:-4]+".txt"
	start = 18
	if "temp" in file:
		end = 54
	else:
		end = 84
	with open(file, 'r') as f_in:
		with open(new_file, 'w') as f_out:
			for i, line in enumerate(f_in):
				if i >= 18 and i <= end:
					f_out.write(line)
