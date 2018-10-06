import os
import glob
import sys
import shutil

cluster_file = sys.argv[1]
with open(cluster_file) as f:
    lines = f.readlines()
lines = [x.strip() for x in lines]

for cluster in lines:
    full_path = cluster.split()[1]
    loop_name = full_path.split('/')[-1]
    folder_name = loop_name[:-4]
    dir_path = '/home/reggiano/loop_grafting/func_loops/dir_cluster_loops/'+folder_name
    os.mkdir(dir_path)
    shutil.copy(full_path, dir_path+'/'+loop_name)
