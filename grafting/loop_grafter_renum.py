#######################################
# old version of loop grafter that    #
# does not relax or design models     #
# only use this if you have a specfic #
# reason                              #
# would not recommend                 #
#######################################


from __future__ import print_function

import argparse
import numpy as np
import pyrosetta
import sys

import pyrosetta.toolbox.numpy_utils as np_utils
from os import path, listdir, makedirs, rename


def configure_insertion_mover(loop, start, stop, overlap=2):
    # Reasons for +1/-1: CCDEndsGraftMover uses an "insert pose into pose" move
    # (written by Steven Lewis) where the "start" residue is the residue before
    # the insertion region, and the "end" residue is the residue after the
    # insertion region; this convention is maintained in this mover.
    # Since we're inserting a whole loop, this puts the required start/end
    # values one residue outside of the loop on either end.
    # Since we define our start/end points as the first/last residues
    # of the loop, we have to subtract/add one to match up with the terminology
    grafter = pyrosetta.rosetta.protocols.grafting.CCDEndsGraftMover(start + 1,
                                                                     stop - 1,
                                                                     loop,
                                                                     overlap,
                                                                     overlap,
                                                                     True)

    # grafter.set_cycles(128)  # defaults to 300

    # scaffold and insert flexibility default to 2, but if there is a larger
    # overlap we should update that here.
    grafter.set_scaffold_flexibility(overlap, overlap)
    grafter.set_insert_flexibility(overlap, overlap)
    grafter.stop_at_closure(False)
    # grafter.copy_pdbinfo(True)
    return grafter


def get_anchor_coordinates_from_pose(p, start, stop):
    """Pack atomic coordinates in a single residue into a numpy matrix.

    Args:
        p (Pose): Pose from coordinates will be extracted
        start (int): the pose-numbered start residue
        stop (int): the pose-numbered stop residue

    Returns:
        Set of coordinates (np.mat) for the atoms in the Pose

    """
    _bb_atoms = ['N', 'CA', 'C', 'O']
    coords = list()
    for resNo in [start, stop]:
        res = p.residue(resNo)

        # only iterate over relevant atoms
        for i in _bb_atoms:
            coords.append([res.xyz(i).x, res.xyz(i).y, res.xyz(i).z])

    return np.mat(coords)


def align_loop_pose_to_anchor_coords(p, coords):
    """Compute and apply the transformation to superpose the residue onto the
    stub.

    Args:
        p (Pose): The Pose instance to manipulate
        coords (np.mat): The coordinates of the atoms in the stub

    Returns:
        The Pose

    """
    moveable_coords = get_anchor_coordinates_from_pose(p, 1, p.size())

    R, t = np_utils.rigid_transform_3D(moveable_coords, coords)
    np_utils.rotate_pose(p, np_utils.numpy_to_rosetta(R))
    np_utils.translate_pose(p, np_utils.numpy_to_rosetta(t.T))

    return p


def insert_loop_into_scaffold(start_pose, loop_fname, loop_range, to_centroid):
    # p = pyrosetta.pose_from_file(scaffold_fname)
    p = start_pose.clone()

    start = loop_range[0]
    stop = loop_range[1]
    
    c = get_anchor_coordinates_from_pose(p, start, stop)
    l = pyrosetta.pose_from_file(loop_fname)
    to_centroid.apply(l)

    looooop = align_loop_pose_to_anchor_coords(l, c)

    # looooop.dump_pdb('transformed_loop.pdb')

    p.delete_residue_range_slow(start, stop)
    smashed = pyrosetta.rosetta.protocols.grafting.insert_pose_into_pose(
        p, looooop, start - 1)

    return smashed   

#added by gabi, renumbering for a scaffold based on given anchors
def renum_pose_and_loop(args, anchors, loop_lengths):
    
    basename,_ = path.splitext(path.basename(args.in_file))
    tmp_file = path.join(args.out_dir,
                                  '{}_tmp.pdb'.format(basename))
                                            
    start_first = anchors[0][0]
    end_first = anchors[0][1]
    current_loop_length_first = end_first - start_first + 1
    new_loop_length_first = loop_lengths[0]
    renum_amt_first = new_loop_length_first - current_loop_length_first
        
    #probably shouldn't be written this way
    if anchors[1] == 0:
            
        renum_bool = False
            
        with open(args.in_file, 'r') as f_in:
            with open(tmp_file, 'w') as f_out:
                for line in f_in:
                    if 'ATOM' in line:
                        if int(line.split()[5]) == end_first and not renum_bool:
                            renum_bool = True
                        if renum_bool:
                            str_old_num = ' '+line.split()[5]+' '
                            str_new_num = ' '+str(int(line.split()[5])+renum_amt_first)+' '
                            if len(str_new_num) - len(str_old_num) == 1:
                                str_new_num = str(int(line.split()[5])+renum_amt_first)+' '
                            new_line = line.replace(str_old_num,str_new_num)
                            f_out.write(new_line)
                        else:
                            if int(line.split()[5]) <= start_first: #this if statement isn't actually necessary
                                f_out.write(line)

    elif len(anchors[1]) == 2:
        
        start_second = anchors[1][0]
        end_second = anchors[1][1]
        current_loop_length_second = end_second - start_second + 1
        new_loop_length_second = loop_lengths[1]
        renum_amt_second = new_loop_length_second - current_loop_length_second
          
        renum_bool_first = False
        renum_bool_second = False
        second_loop_encountered = 2
           
        with open(args.in_file, 'r') as f_in:
            with open(tmp_file, 'w') as f_out:
                for line in f_in:
                    if 'ATOM' in line:
                        if int(line.split()[5]) == end_first and not renum_bool_first:
                            renum_bool_first = True
                        if int(line.split()[5]) == end_second and not  renum_bool_second:
                            renum_bool_second = True
                        if renum_bool_first and renum_bool_second:
                            str_old_num = ' '+line.split()[5]+' '
                            str_new_num = ' '+str(int(line.split()[5])+renum_amt_first+renum_amt_second)+' '
                            if len(str_new_num) - len(str_old_num) == 1:
                                str_new_num = str(int(line.split()[5])+renum_amt_first+renum_amt_second)+' '
                            new_line = line.replace(str_old_num,str_new_num)
                            f_out.write(new_line)
                        elif renum_bool_first:
                            if second_loop_encountered == 2:
                                second_loop_encountered = 1
                            str_old_num = ' '+line.split()[5]+' '
                            str_new_num = ' '+str(int(line.split()[5])+renum_amt_first)+' '
                            if len(str_new_num) - len(str_old_num) == 1:
                                str_new_num = str(int(line.split()[5])+renum_amt_first)+' '
                            new_line = line.replace(str_old_num,str_new_num)
                            f_out.write(new_line)
                        elif renum_bool_second:
                            if second_loop_encountered == 2:
                                second_loop_encountered = 0
                            str_old_num = ' '+line.split()[5]+' '
                            str_new_num = ' '+str(int(line.split()[5])+renum_amt_second)+' '
                            if len(str_new_num) - len(str_old_num) == 1:
                                str_new_num = str(int(line.split()[5])+renum_amt_second)+' '
                            new_line = line.replace(str_old_num,str_new_num)
                            f_out.write(new_line)
                        else:
                            f_out.write(line)
            
        #loop range of second loop also has to be renumbered
        if second_loop_encountered == 0: #first anchor pair was seen 2nd
            
            new_start = anchors[0][0] + renum_amt_second
            new_end = anchors[0][1] + renum_amt_second
            for i, directory in enumerate(args.loop_dir):
                if i == 0:
                    components = directory.split('_')
                    #needs editing to get name
                    no_range = components[0]
                    new_dir = no_range+'_'+str(new_start)+'_'+str(new_end)
                    rename(directory, new_dir)
                    args.loop_dir[0] = new_dir
                    
        elif second_loop_encountered == 1: #second anchor pair was seen 2nd

            new_start = anchors[1][0] + renum_amt_first
            new_end = anchors[1][1] + renum_amt_first
            for i, directory in enumerate(args.loop_dir):
                if i == 1:
                    components = directory.split('_')
                    #needs editing to get old range
                    no_range = components[0]
                    new_dir = no_range+'_'+str(new_start)+'_'+str(new_end)
                    rename(directory, new_dir)
                    args.loop_dir[1] = new_dir
    
    p = pyrosetta.pose_from_file(tmp_file)
    return p
    
def main(argv):
    parser = argparse.ArgumentParser(description='Program')
    parser.add_argument('-i', '--in_file', action='store', type=str,
                        help='input PDB scaffold file')
    parser.add_argument('-d', '--loop_dir', action='store', type=str,
                        required=True, nargs='+',
                        help='directories containing loops to graft. Provide'
                        ' directories with C-terminal-most loops first.')
    parser.add_argument('-o', '--out_dir', action='store', type=str,
                        help='directory to write output')
    #added by gabi
    parser.add_argument('-a1', '--anchor_res_pair1', nargs='+', type=int,
                        help='first pair of anchor residues')
    parser.add_argument('-a2', '--anchor_res_pair2', nargs='+', type=int,
                        help='second pair of anchor residues', default=0)
    args = parser.parse_args()
    
#    print(list_of_loop_files[0])
    
    #added by gabi
    anchor_list = [args.anchor_res_pair1, args.anchor_res_pair2]
    if not path.exists(args.out_dir):
        makedirs(args.out_dir)

    # user-specified scorefxn
    pyrosetta.init()  # extra_options='-beta_nov15')

    # figure out the pdb-numbered residue ranges for the loops from the
    # directory name. This is hacky, but whateva
    def _get_loop_ranges(dir_name, split_char='_'):
        components = dir_name.split(split_char)
        res_range = list()
        for c in components:
            try:
                res_range.append(int(c))
            except ValueError:
                pass
        return res_range

    def _get_loop_size(dir_name):
        dir_path = path.join(dir_name)
        loop_file = next(path.join(dir_path,f) for f in listdir(dir_path) if path.isfile(path.join(dir_path, f)))
        loop_pose = pyrosetta.pose_from_file(loop_file)
        loop_size = loop_pose.total_residue()
        return loop_size

    loop = list()
    start_pose = pyrosetta.pose_from_file(args.in_file)
    #added by gabi
    loops_length = list()
    for directory in args.loop_dir:
        loop_res = _get_loop_size(directory)
        loops_length.append(loop_res)
    edited_pose = renum_pose_and_loop(args, anchor_list, loops_length)
    
    #pose_anchor_res = get_anchor_for_renum(args,len(list_of_loop_files))
    #moved down here by gabi
    list_of_loop_files = [[path.join(ld, l) for l in
                          listdir(ld)] for ld in args.loop_dir]
    to_centroid = pyrosetta.rosetta.protocols.simple_moves.SwitchResidueTypeSetMover('centroid')
    to_centroid.apply(edited_pose)

    for size in loops_length:
        
        loop_range = _get_loop_ranges(directory)
        loop.append((edited_pose.pdb_info().pdb2pose('A', loop_range[0]),
                     edited_pose.pdb_info().pdb2pose('A', loop_range[1])))
        
        assert(loop[-1][0] < loop[-1][1] and loop[-1][0] > 0)

    # temporary -- I'd rather have some recursive type of thing that doesn't
    # need to know how many loops we're talking about here
    assert(len(list_of_loop_files) ==  len(loop) )
    print(len(list_of_loop_files))
    sf = pyrosetta.get_score_function(False)  # centroid score function
    score_manager = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
    score_term = score_manager.score_type_from_name("atom_pair_constraint")
    sf.set_weight(score_term, 1.0)

    #added by gabi
    score_term = score_manager.score_type_from_name("vdw")
    sf.set_weight(score_term, 1.0)
    score_term = score_manager.score_type_from_name("hbond_sr_bb")
    sf.set_weight(score_term, 1.0)
    score_term = score_manager.score_type_from_name("hbond_lr_bb")
    sf.set_weight(score_term, 0.5)
    score_term = score_manager.score_type_from_name("dihedral_constraint")
    sf.set_weight(score_term, 10.0)
    score_term = score_manager.score_type_from_name("rama")
    sf.set_weight(score_term, 0.1)
    score_term = score_manager.score_type_from_name("omega")
    sf.set_weight(score_term, 1.0)
    score_term = score_manager.score_type_from_name("rg")
    sf.set_weight(score_term, 1.0)
    score_term = score_manager.score_type_from_name("linear_chainbreak")
    sf.set_weight(score_term, 2.778)
    score_term = score_manager.score_type_from_name("cart_bonded")
    sf.set_weight(score_term, 0.5)

    print(sf)
    # to avoid bookkeeping crap, we're going to smash in the downstream loop
    # first.
    if len(list_of_loop_files) == 1:
        print("One loop!")
        for i, loop_fname in enumerate(list_of_loop_files[0]):
            p = insert_loop_into_scaffold(edited_pose, loop_fname, loop[0], to_centroid)
            
            #try to filter based on the cartbonded score, misses one or two usually
            if sf.score_by_scoretype(p,score_term,True) < 180.:
                scaffold_base, _ = path.splitext(path.basename(args.in_file))
                l1_basename, _ = path.splitext(path.basename(loop_fname))
                fname = path.join(args.out_dir,
                                  '{}_{}.pdb'.format(scaffold_base, l1_basename))
                p.dump_pdb(fname)

    elif len(list_of_loop_files) == 2:
        for i, loop_fname in enumerate(list_of_loop_files[0]):
            p = insert_loop_into_scaffold(edited_pose, loop_fname, loop[0],
                                      to_centroid)
            for i, loop2_fname in enumerate(list_of_loop_files[1]):
                p2 = insert_loop_into_scaffold(p, loop2_fname, loop[1], to_centroid)
                
                #may need to work on the best filter score
                if sf.score_by_scoretype(p,score_term,True) < 180.:
                    scaffold_base, _ = path.splitext(path.basename(args.in_file))
                    l1_basename, _ = path.splitext(path.basename(loop_fname))
                    l2_basename, _ = path.splitext(path.basename(loop2_fname))
                    fname = path.join(args.out_dir,
                                  '{}_{}_{}.pdb'.format(scaffold_base,
                                                        l1_basename,
                                                        l2_basename))
                    p2.dump_pdb(fname)


if __name__ == '__main__':
    main(sys.argv)
