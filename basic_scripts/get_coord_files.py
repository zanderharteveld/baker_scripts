################################################################
# this script creates a coordinate file for the specified loop #
# all other residues will be constrained so that only the loop #
# is relaxed                                                   #
################################################################

import os
import glob

for pdb in glob.glob("*.pdb"):
    coord_file = pdb.strip('.pdb') + '_coords.cst'
    c_type = 'CoordinateConstraint'
    a_name = 'CA'
    loop_to_graft = [95,98]
    with open(pdb, 'r') as f_in:
        with open(coord_file, 'w') as f_out:
            refatom = 0
            for line in f_in:
                if 'CA' in line:
                    line = line.strip().split()
                    if refatom == 0:
                        refatom = int(line[5])
                    if int(line[5]) < loop_to_graft[0] or int(line[5]) > loop_to_graft[1]:
                        resnum = int(line[5])
                        coord1 = float(line[6])
                        coord2 = float(line[7])
                        coord3 = float(line[8])
                        new_line = '%s %s %s %s %s %s %s %s HARMONIC 0.0 0.2 \n' %(c_type,a_name,resnum,a_name,refatom,coord1,coord2,coord3)
                        f_out.write(new_line)
