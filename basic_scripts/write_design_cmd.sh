#!/bin/bash

#########################################################################################
# this script is to write rosetta scripts command line                                  #
# when a variable is involved for each input pdb                                        #
# call this script with:                                                                #
# for i in out/segment_profile/*pdb; do sh write_design_cmd.sh $i; done  > design_tasks #
#########################################################################################

i=$1

j=${i/.pdb/}
j+=".resfile"
echo "~/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease -parser:protocol design_loop.xml -parser:script_vars res_file=${j} -out:path:all out/design_hairpin/trial4 -s ${i} -database /software/rosetta/main/database/ -overwrite -nstruct 5"
