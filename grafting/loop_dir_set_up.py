import os
import shutil
#######################################################
# this script will put one loop in it's own directory #
# inside of an output directory, use if want to do    #
# fake_parallel_graft                                 #
# inputs are the directory containing all loop pdbs   #
# and the name of the output directory                #
#######################################################


import glob
import sys
import argparse
import fileinput

input_dir = sys.argv[1]
output_dir = sys.argv[2]

os.chdir(input_dir)
input_path = os.getcwd()
pdbs = sorted(glob.glob('*.pdb'))
os.chdir('../')
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
os.chdir(output_dir)
output_path = os.getcwd()
for i in range(0,len(pdbs)-2):
    loop1 = pdbs[i]
    folder_name = loop1[:-4]
    os.mkdir(folder_name)
    os.mkdir(folder_name+'/single')
    os.mkdir(folder_name+'/combo')
    single_dir_fullpath = output_path+'/'+folder_name+'/single'
    loop1_fullpath = input_path+'/'+loop1
    shutil.copyfile(loop1_fullpath,single_dir_fullpath+'/'+loop1)
    j = i+1
    while j < len(pdbs):
        loop2 = pdbs[j]
      #  if j+1 >= len(pdbs):
      #      loop3 = 0
      #  else:
      #      loop3 = pdbs[j+1]
        # checks that loop names are same
        # will make sure that no loops of different sizes are put in the same directory
        # would mess up the loop_grafter
      #  if loop3 != 0 and loop2[0:3] == loop3[0:3]:
      #      combo_dir_name = 'combine_'loop2[:-4]+'_'+loop3[:-4]
      #      os.mkdir(folder_name+'/combo/'+combo_dir_name)
      #      combo_dir_fullpath = output_path+'/'+folder_name+'/combo/'+combo_dir_name
      #      loop2_fullpath = input_path+'/'+loop2
      #      loop3_fullpath = input_path+'/'+loop3
      #      shutil.copyfile(loop2_fullpath,combo_dir_fullpath+'/'+loop2)
      #      shutil.copyfile(loop3_fullpath,combo_dir_fullpath+'/'+loop3)
      #      j += 2
      #  else:
        combo_dir_name = loop2[:-4]
        os.mkdir(folder_name+'/combo/'+combo_dir_name)
        combo_dir_fullpath = output_path+'/'+folder_name+'/combo/'+combo_dir_name
        loop2_fullpath = input_path+'/'+loop2
        shutil.copyfile(loop2_fullpath,combo_dir_fullpath+'/'+loop2)
        j += 1



    

