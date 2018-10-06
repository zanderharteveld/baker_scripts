################################################
# this script uses loop grafter to             #
# design new loops for shape complementarity   #
# there are also options for no SC             #
# renumbering of loops and pdbs not neccessary #
# must give C-terminal most loop first         #
# need flags file and score xml in cwd         #
# both can be found in:                        #
# /home/reggiano/scripts/xml                   #
################################################

from __future__ import print_function

import argparse
import numpy as np
import pyrosetta
import sys
import shutil

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
def renum_pose_and_loop(args, anchors, loop_ranges):
    
    basename,_ = path.splitext(path.basename(args.in_file))
    tmp_file = path.join(args.out_dir,
                                  '{}_tmp.pdb'.format(basename))
                                            
    start_first = anchors[0][0]
    end_first = anchors[0][1]
    current_loop_length_first = end_first - start_first
    new_loop_length_first = loop_ranges[0][1] - loop_ranges[0][0]
    renum_amt_first = new_loop_length_first - current_loop_length_first

    #probably shouldn't be written this way
    if anchors[1] == 0:
        
        ter = False   
        renum_bool = False
            
        with open(args.in_file, 'r') as f_in:
            with open(tmp_file, 'w') as f_out:
                for line in f_in:
                    if args.shape_complement:
                        if 'HETATM' in line:
                            if not ter:
                                f_out.write('TER\n')
                                ter = True
                            str_old_num = ' '+line.split()[5]+' '
                            str_new_num = ' '+str(int(line.split()[5])+renum_amt_first)+' '
                            new_line = line.replace(str_old_num,str_new_num)
                            f_out.write(new_line)
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
        current_loop_length_second = end_second - start_second
        new_loop_length_second = loop_ranges[1][1] - loop_ranges[1][0]
        renum_amt_second = new_loop_length_second - current_loop_length_second
          
        renum_bool_first = False
        renum_bool_second = False
        ter = False
           
        with open(args.in_file, 'r') as f_in:
            with open(tmp_file, 'w') as f_out:
                for line in f_in:
                    if args.shape_complement:
                        if 'HETATM' in line:
                            if not ter:
                                f_out.write('TER\n')
                                ter = True
                            str_old_num = ' '+line.split()[5]+' '
                            str_new_num = ' '+str(int(line.split()[5])+renum_amt_first+renum_amt_second)+' '
                            new_line = line.replace(str_old_num,str_new_num)
                            f_out.write(new_line)
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
                            str_old_num = ' '+line.split()[5]+' '
                            str_new_num = ' '+str(int(line.split()[5])+renum_amt_first)+' '
                            if len(str_new_num) - len(str_old_num) == 1:
                                str_new_num = str(int(line.split()[5])+renum_amt_first)+' '
                            new_line = line.replace(str_old_num,str_new_num)
                            f_out.write(new_line)
                        elif renum_bool_second:
                            str_old_num = ' '+line.split()[5]+' '
                            str_new_num = ' '+str(int(line.split()[5])+renum_amt_second)+' '
                            if len(str_new_num) - len(str_old_num) == 1:
                                str_new_num = str(int(line.split()[5])+renum_amt_second)+' '
                            new_line = line.replace(str_old_num,str_new_num)
                            f_out.write(new_line)
                        else:
                            f_out.write(line)
    
    if args.shape_complement:
        p = pyrosetta.rosetta.core.pose.Pose()
        ligand_params = pyrosetta.Vector1([args.ligand])
        res_set = p.conformation().modifiable_residue_type_set_for_conf()
        res_set.read_files_for_base_residue_types(ligand_params)
        p.conformation().reset_residue_type_set_for_conf(res_set)
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
                        help='first pair of anchor residues, where the new loop '
                            'should be grafted')
    parser.add_argument('-a2', '--anchor_res_pair2', nargs='+', type=int,
                        help='second pair of anchor residues', default=0)
    parser.add_argument('-sc', '--shape_complement', type=bool, default=False,
                        help='do you want to optimize loops for shape complementarity'
                            ' of a ligand?')
    parser.add_argument('-l', '--ligand', default = 'LIG.fa.params', type=str,
                        help='ligand params file')
    args = parser.parse_args()
    
#    print(list_of_loop_files[0])
    
    if not path.exists(args.out_dir):
        makedirs(args.out_dir)
    
    #for ligands
    if not args.shape_complement:
        with open('./flags', 'r') as f:
            lines = f.readlines()
            pyrosetta.init(' '.join(lines))
    else:
        option = '-extra_res_fa '+args.ligand
        with open('./flags', 'r') as f:
            lines = f.readlines()
            pyrosetta.init(' '.join(lines)+' '+option)
    
    #get length of loop based on first file in loop dir, all loops must be same length
    def _get_loop_length(dir_name):
        dir_path = path.join(dir_name)
        loop_file = next(path.join(dir_path,f) for f in listdir(dir_path) if path.isfile(path.join(dir_path, f)))
        loop_pose = pyrosetta.pose_from_file(loop_file)
        loop_size = loop_pose.total_residue()
        return loop_size
    
    loop = list()
    loop_anchors = list()
    loop_lengths = list()
    
    #use loop length and orig anchors to get new renumbered anchors
    #assumes that C-terminal most loop is given first
    orig_anchors = [args.anchor_res_pair1, args.anchor_res_pair2]
    for directory in args.loop_dir:
        loop_res = _get_loop_length(directory)
        loop_lengths.append(loop_res)
    
    if len(loop_lengths) == 1:
        orig_length = orig_anchors[0][1]-orig_anchors[0][0]+1
        renum_amt = loop_lengths[0] - orig_length
        loop_anchors.append((orig_anchors[0][0], orig_anchors[0][1]+renum_amt))
        loop_anchors.append(0)
    else:
        orig_length_C_term = orig_anchors[0][1]-orig_anchors[0][0]+1
        orig_length_N_term = orig_anchors[1][1]-orig_anchors[1][0]+1
        renum_amt_C_term = loop_lengths[0] - orig_length_C_term
        renum_amt_N_term = loop_lengths[1] - orig_length_N_term
        loop_anchors.append((orig_anchors[0][0]+renum_amt_N_term,
                                orig_anchors[0][1]+renum_amt_N_term+renum_amt_C_term))
        loop_anchors.append((orig_anchors[1][0], orig_anchors[1][1]+renum_amt_N_term))
         
    edited_pose = renum_pose_and_loop(args, orig_anchors, loop_anchors)
    
    #for residue selectors
    index_str = ''
    if loop_anchors[1] == 0:
            for res in range(loop_anchors[0][0], loop_anchors[0][1]+1):
                if res != loop_anchors[0][1]:
                    index_str+=str(res)+','
                else:
                    index_str+=str(res)
    else:
        for pair in loop_anchors:
            for res in range(pair[0], pair[1]+1):
                if res != loop_anchors[1][1]:
                    index_str+=str(res)+','
                else:
                    index_str+=str(res)
    
    
    #moved down here by gabi
    list_of_loop_files = [[path.join(ld, l) for l in
                          listdir(ld)] for ld in args.loop_dir]
    to_centroid = pyrosetta.rosetta.protocols.simple_moves.SwitchResidueTypeSetMover('centroid')
    to_centroid.apply(edited_pose)

    for loop_range in loop_anchors:
        loop.append((edited_pose.pdb_info().pdb2pose('A', loop_range[0]),
                     edited_pose.pdb_info().pdb2pose('A', loop_range[1])))
        assert(loop[-1][0] < loop[-1][1] and loop[-1][0] > 0)

    # temporary -- I'd rather have some recursive type of thing that doesn't
    # need to know how many loops we're talking about here
    assert(len(list_of_loop_files) ==  len(loop) )
    print(len(list_of_loop_files))
    sf = pyrosetta.get_score_function(False)  # centroid score function
    fa_sf = pyrosetta.create_score_function('ref2015_cart')   # for filtering after relax/mcdesign
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
    score_term = score_manager.score_type_from_name("cart_bonded")
    sf.set_weight(score_term, 0.5)
    score_term = score_manager.score_type_from_name("linear_chainbreak")
    sf.set_weight(score_term, 2.778)

    print(sf)
    
    #set up generalized min_mover and design stuff
    switch_mover = pyrosetta.rosetta.protocols.simple_moves.SwitchResidueTypeSetMover('full_atom')

    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_bb(False)
    mm.set_chi(False)

    tf_fd = pyrosetta.rosetta.core.pack.task.TaskFactory()
    tf_fd.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf_fd.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf_fd.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

    tf_mp = pyrosetta.rosetta.core.pack.task.TaskFactory()
    tf_mp.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf_mp.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf_mp.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

    # to avoid bookkeeping crap, we're going to smash in the downstream loop
    # first.
    if len(list_of_loop_files) == 1:
        print("One loop!")
        for i, loop_fname in enumerate(list_of_loop_files[0]):
            p = insert_loop_into_scaffold(edited_pose, loop_fname, loop[0], to_centroid)
            
            #p.dump_pdb("centroid.pdb")
            switch_mover.apply(p)
            #recover old rotamers, not working, nres == init_pose->size() failed
            #save_sc_mover.apply(p)

            #finish set up for minmover
            for res in range(loop_anchors[0][0]-3, loop_anchors[0][1]+4):
                mm.set_bb(res,True)
                mm.set_chi(res,True)

            chainB = pyrosetta.rosetta.core.select.residue_selector.ChainSelector("B")
            tf_mp.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
                pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), chainB))
            tf_mp.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
                pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT(), chainB))
            tf_mp.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
                pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), chainB, True))
            task = tf_mp.create_task_and_apply_taskoperations(p)
            pack_mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(fa_sf,task)
            min_mover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
            min_mover.cartesian(True)
            min_mover.set_movemap(mm)
            min_mover.tolerance(1E-4)
            min_mover.score_function(fa_sf)
            pack_mover.apply(p)
            min_mover.apply(p)

            cb_filter = pyrosetta.rosetta.protocols.simple_filters.ChainBreak()
            passed = cb_filter.apply(p)

            #p.dump_pdb("quick_min.pdb")

            if passed:
                if args.shape_complement:
                    #do set up for task operations
                    loop1 = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(index_str)
                    not_loop = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(loop1)
                    tf_fd.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
                        pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT(), not_loop))
                    no_M_C = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
                    no_M_C.aas_to_keep("RHKDESTNQGPAVILFYW")
                    tf_fd.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(no_M_C, both_loops))
                    tf_fd.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
                        pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), not_loop))
                    fast_design = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign(fa_sf,1)
                    fast_design.set_task_factory(tf_fd)
                    fast_design.cartesian(True)
                    fast_design.min_type("lbfgs_armijo_nonmonotone")
                    fast_design.set_movemap(mm)
                    sc_filter = pyrosetta.rosetta.protocols.simple_filters.ShapeComplementarityFilter(0.,0.,1,0,0)
                    tot_score = score_manager.score_type_from_name("total_score")
                    tot_score_filter = pyrosetta.rosetta.protocols.simple_filters.ScoreTypeFilter(fa_sf,tot_score,0)
                    #gmc arguments: maxtrials, max accepted trials, task scaling, mover, temperature, sample_type, drift
                    gmc = pyrosetta.rosetta.protocols.simple_moves.GenericMonteCarloMover(5,1,0,fast_design,1.,"low",False)
                    #add_filter arguments: filter, adaptive, temp, sample_type
                    gmc.add_filter(sc_filter, False, 1., "high")
                    gmc.add_filter(tot_score_filter, False, 1., "low")
                    gmc.apply(p2)

                scaffold_base, _ = path.splitext(path.basename(args.in_file))
                l1_basename, _ = path.splitext(path.basename(loop_fname))
                fname = path.join(args.out_dir,
                                  '{}_{}_{}_{}.pdb'.format(scaffold_base, l1_basename, 
						              loop_anchors[0][0], loop_anchors[0][1]))
                scoring_xml = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_file("./score.xml")
                protocol = scoring_xml.get_mover("ParsedProtocol")
                protocol.apply(p)
                p.dump_pdb(fname)

    elif len(list_of_loop_files) == 2:
        for i, loop_fname in enumerate(list_of_loop_files[0]):
            p = insert_loop_into_scaffold(edited_pose, loop_fname, loop[0],
                                      to_centroid)
            for i, loop2_fname in enumerate(list_of_loop_files[1]):
                p2 = insert_loop_into_scaffold(p, loop2_fname, loop[1], to_centroid)               
                
                #p2.dump_pdb("centroid.pdb")
                switch_mover.apply(p2)
                #revover old rotamers
                #save_sc_mover.apply(p2)

                #finish set up of movemap
                for pair in loop_anchors:
                    for res in range(pair[0]-3,pair[1]+4): #add in 3 residues to minimize
                        mm.set_bb(res,True)
                        mm.set_chi(res,True)

                #reduce clashes before do fast design, filter out ones that are no good afterwards
                chainB = pyrosetta.rosetta.core.select.residue_selector.ChainSelector("B")
                tf_mp.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
                    pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), chainB))
                tf_mp.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
                    pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT(), chainB))
                tf_mp.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
                    pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), chainB, True))
                task = tf_mp.create_task_and_apply_taskoperations(p2)
                pack_mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(fa_sf,task)
                min_mover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
                min_mover.cartesian(True)
                min_mover.set_movemap(mm)
                min_mover.tolerance(1E-4)
                min_mover.score_function(fa_sf)
                pack_mover.apply(p2)
                min_mover.apply(p2)

                #p2.dump_pdb("quick_min.pdb")

                cb_filter = pyrosetta.rosetta.protocols.simple_filters.ChainBreak()
                passed = cb_filter.apply(p2)       

                #pack and minize before design, take out shitty ones
                if passed:
                    if args.shape_complement:
                        #do set up for task operations
                        both_loops = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(index_str)
                        not_loop = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(both_loops)
                        tf_fd.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
                            pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT(), not_loop))
                        no_M_C = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
                        no_M_C.aas_to_keep("RHKDESTNQGPAVILFYW")
                        tf_fd.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(no_M_C, both_loops))
                        tf_fd.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
                            pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), not_loop))
                        fast_design = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign(fa_sf,1)
                        fast_design.set_task_factory(tf_fd)
                        fast_design.cartesian(True)
                        fast_design.min_type("lbfgs_armijo_nonmonotone")
                        fast_design.set_movemap(mm)
                        sc_filter = pyrosetta.rosetta.protocols.simple_filters.ShapeComplementarityFilter(0.,0.,1,0,0)
                        tot_score = score_manager.score_type_from_name("total_score")
                        tot_score_filter = pyrosetta.rosetta.protocols.simple_filters.ScoreTypeFilter(fa_sf,tot_score,0)
                        #gmc arguments: maxtrials, max accepted trials, task scaling, mover, temperature, sample_type, drift
                        gmc = pyrosetta.rosetta.protocols.simple_moves.GenericMonteCarloMover(5,1,0,fast_design,1.,"low",False)
                        #add_filter arguments: filter, adaptive, temp, sample_type
                        gmc.add_filter(sc_filter, False, 1., "high")
                        gmc.add_filter(tot_score_filter, False, 1., "low")
                        gmc.apply(p2) 

                    scaffold_base, _ = path.splitext(path.basename(args.in_file))
                    l1_basename, _ = path.splitext(path.basename(loop_fname))
                    l2_basename, _ = path.splitext(path.basename(loop2_fname))
                    fname = path.join(args.out_dir,
                                  '{}_{}_{}_{}_{}_{}_{}.pdb'.format(scaffold_base,
                                                        l1_basename, loop_anchors[0][0], loop_anchors[0][1],
                                                        l2_basename, loop_anchors[1][0], loop_anchors[1][1]))
                    scoring_xml = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_file("./score.xml")
                    protocol = scoring_xml.get_mover("ParsedProtocol")
                    protocol.apply(p2)
                    p2.dump_pdb(fname)

if __name__ == '__main__':
    main(sys.argv)
