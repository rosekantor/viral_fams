#!/usr/bin/python3

import sys
import os
import argparse
import re
import subprocess
from Bio import Phylo
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import Normalize
from collections import defaultdict

# software run
MUSCLE_bin = "muscle"
GBLOCKS_bin = "Gblocks"
FASTTREE_bin = "FastTree"


# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="make trees from sequences for fam hits")

    parser.add_argument("-d", "--datadir", help="")
    parser.add_argument("-r", "--ref-seq-dir", help="")
    parser.add_argument("-c", "--coverage", default=0.5, type=float, help="fraction of overlap region required covered [0.0,1.0] (def: 0.5)")
    parser.add_argument("-b", "--bitscore-map", help="")
    parser.add_argument("-o", "--outdir", help="")
    
    args = parser.parse_args()
    args_pass = True

    if args.datadir is None:
        print ("must specify --{}\n".format('datadir'))
        args_pass = False
    elif not os.path.exists(args.datadir) or \
         not os.path.isdir(args.datadir):
        print ("--{} {} must exist\n".format('datadir', args.datadir))
        args_pass = False

    if args.ref_seq_dir is None:
        print ("must specify --{}\n".format('ref-seq-dir'))
        args_pass = False
    elif not os.path.exists(args.ref_seq_dir) or \
         not os.path.isdir(args.ref_seq_dir):
        print ("--{} {} must exist\n".format('ref-seq-dir', args.ref_seq_dir))
        args_pass = False

    if args.outdir is None:
        print ("must specify --{}\n".format('outdir'))
        args_pass = False
    elif not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    elif not os.path.isdir(args.outdir):
        print ("--{} {} is not a dir\n".format('outdir', args.outdir))
    
    if args.bitscore_map:
        if not os.path.exists(args.bitscore_map):
            print ("--{} {} must exist\n".format('bitscore-map', args.bitscore_map))
        
    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


# get_fam_seq_path ()
#
def get_fam_seq_path (data_dir, faa_file):
    """
    fam_dir = os.path.join (data_dir, fam_id)
    fam_seq_file = fam_id+'-sliced.faa'
    #fam_seq_file = fam_id+'-no_slice_label-sliced.faa'
    #fam_seq_file = fam_id+'-sp_label-sliced.faa'
    fam_seq_path = os.path.join (fam_dir, fam_seq_file)
    """
    fam_seq_path = os.path.join(data_dir, os.path.basename(faa_file))
    return fam_seq_path


# get_fam_seq_nr_path ()
#
def get_fam_seq_nr_path (data_dir, faa_file):
    """
    fam_dir = os.path.join (data_dir, fam_id)
    fam_seq_file = fam_id+'-sliced.faa'
    #fam_seq_file = fam_id+'-no_slice_label-sliced.faa'
    #fam_seq_file = fam_id+'-sp_label-sliced.faa'
    fam_seq_path = os.path.join (fam_dir, fam_seq_file)
    """
    fam_nr_seq_path = os.path.join(data_dir, os.path.basename(faa_file).replace('.faa','-nr.faa'))
    return fam_nr_seq_path


# get_muscle_output_path ()
#
def get_muscle_output_path (data_dir, faa_file):
    """
    fam_dir = os.path.join (data_dir, fam_id)
    muscle_file = fam_id+'-sliced-muscle.afa'
    #muscle_file = fam_id+'-no_slice_label-sliced-muscle.afa'
    #muscle_file = fam_id+'-sp_label-sliced-muscle.afa'
    muscle_path = os.path.join (fam_dir, fam_seq_file)
    """
    muscle_path = os.path.join (data_dir, os.path.basename(faa_file).replace('.faa','-nr-muscle.afa'))
    return muscle_path

# get_filtered_list_output_path ()
#
def get_filtered_list_output_path (data_dir, faa_file):
    """
    fam_dir = os.path.join (data_dir, fam_id)
    muscle_file = fam_id+'-sliced-muscle.afa'
    #muscle_file = fam_id+'-no_slice_label-sliced-muscle.afa'
    #muscle_file = fam_id+'-sp_label-sliced-muscle.afa'
    muscle_path = os.path.join (fam_dir, fam_seq_file)
    """
    filtered_list_output_path = os.path.join (data_dir, os.path.basename(faa_file).replace('.faa','-nr-muscle-filtered-out.txt'))
    return filtered_list_output_path
                                   
# get_muscle_filtered_output_path ()
#
def get_muscle_filtered_output_path (data_dir, faa_file):
    """
    fam_dir = os.path.join (data_dir, fam_id)
    muscle_file = fam_id+'-sliced-muscle.afa'
    #muscle_file = fam_id+'-no_slice_label-sliced-muscle.afa'
    #muscle_file = fam_id+'-sp_label-sliced-muscle.afa'
    muscle_path = os.path.join (fam_dir, fam_seq_file)
    """
    muscle_filtered_path = os.path.join (data_dir, os.path.basename(faa_file).replace('.faa','-nr-muscle-filt.afa'))
    return muscle_filtered_path


# get_gblocks_output_path ()
#
def get_gblocks_output_path (data_dir, faa_file):
    """
    fam_dir = os.path.join (data_dir, fam_id)
    gblocks_file = fam_id+'-sliced-muscle-gblocks.afa'
    #gblocks_file = fam_id+'-no_slice_label-sliced-muscle-gblocks.afa'
    #gblocks_file = fam_id+'-sp_label-sliced-muscle-gblocks.afa'
    gblocks_path = os.path.join (fam_dir, fam_seq_file)
    """
    gblocks_path = os.path.join (data_dir, os.path.basename(faa_file).replace('.faa','-nr-muscle-filt-gblocks.afa'))
    return gblocks_path

# get_modified_msa_output_path ()
#
def get_modified_msa_output_path (data_dir, faa_file):
    """
    fam_dir = os.path.join (data_dir, fam_id)
    mod_gblocks_file = fam_id+'-sliced-muscle-gblocks-modified.afa'
    #mod_gblocks_file = fam_id+'-no_slice_label-sliced-muscle-gblocks-modifield.afa'
    #mod_gblocks_file = fam_id+'-sp_label-sliced-muscle-gblocks-modifield.afa'
    mod_gblocks_path = os.path.join (fam_dir, fam_seq_file)
    """
    mod_gblocks_path = os.path.join (data_dir, os.path.basename(faa_file).replace('.faa','-nr-muscle-filt-gblocks-modified.afa'))
    return mod_gblocks_path


# get_tree_output_path ()
#
def get_tree_output_path (data_dir, faa_file):
    """
    fam_dir = os.path.join (data_dir, fam_id)
    newick_tree_file = fam_id+'-sliced-muscle-gblocks-fasttree.newick'
    #newick_tree_file = fam_id+'-no_slice_label-sliced-muscle-gblocks-fasttree.newick'
    #newick_tree_file = fam_id+'-sp_label-sliced-muscle-gblocks-fasttree.newick'
    newick_tree_path = os.path.join (fam_dir, newick_tree_file)
    """
    newick_tree_path = os.path.join (data_dir, os.path.basename(faa_file).replace('.faa','-nr-muscle-filt-gblocks-fasttree.newick'))
    return newick_tree_path


# get_id_map_path ()
#
def get_id_map_path (data_dir, faa_file):
    """
    fam_dir = os.path.join (data_dir, fam_id)
    id_map_file = fam_id+'-sliced-muscle-gblocks-fasttree.id_map'
    #id_map_file = fam_id+'-no_slice_label-sliced-muscle-gblocks-fasttree.id_map'
    #id_map_file = fam_id+'-sp_label-sliced-muscle-gblocks-fasttree.id_map'
    id_map_path = os.path.join (fam_dir, id_map_file)
    """
    id_map_path = os.path.join (data_dir, os.path.basename(faa_file).replace('.faa','-nr-muscle-filt-gblocks-fasttree.id_map'))
    return id_map_path

# get_phylo_vis_output_path(data_dir, faa_file)
#
def get_phylo_vis_output_path(data_dir, faa_file):
    phylo_vis_output_path = os.path.join (data_dir, os.path.basename(faa_file).replace('.faa','-nr-muscle-filt-gblocks-fasttree-vis.png'))
    return phylo_vis_output_path


# read_fam_ids()
#
def read_fam_ids (fam_id_listfile):
    fam_ids = []
    
    with open (fam_id_listfile, 'r') as fid_h:
        for fid_line in fid_h:
            fid_line = fid_line.strip()
            if not fid_line or fid_line.startswith('#'):
                continue
            rec = fid_line.split("\t")
            fam_id = rec[0]
            fam_ids.append(fam_id)
                
    return fam_ids


# make_short_id_msa ()
#
def make_short_id_msa (muscle_msa_path):
    short_id_map = dict()
    short_id_msa_path = muscle_msa_path+'-short_ids'

    outbuf = []
    
    index = 0
    with open (muscle_msa_path,'r') as msa_orig:
        for m_line in msa_orig:
            m_line = m_line.strip()
            if m_line.startswith ('>'):
                header = m_line
                index += 1 
                short_id = 'S_'+str(index)  
                short_id_map[short_id] = header
                outbuf.append('>'+short_id)
            else:
                outbuf.append(m_line)

    with open (short_id_msa_path,'w') as msa_new:
        msa_new.write("\n".join(outbuf)+"\n")
        
    return (short_id_msa_path, short_id_map)


# fix_gblocks_output ()
#
def fix_gblocks_output (short_id_msa_gblocks_out_path,
                        gblocks_output_path,
                        short_id_map):

    outbuf = []
    seq = ''
    
    with open (short_id_msa_gblocks_out_path,'r') as short_msa:
        for msa_line in short_msa:
            msa_line = msa_line.strip()
            if msa_line.startswith ('>'):
                if seq:
                    outbuf.append(seq)
                    seq = ''
                short_id = msa_line.split()[0].replace('>','')
                outbuf.append(short_id_map[short_id])
            else:
                seq += msa_line.replace(' ','')
        if seq:
            outbuf.append(seq)

    with open (gblocks_output_path,'w') as long_msa:
        long_msa.write("\n".join(outbuf)+"\n")

    return 

# get_modified_id(id)
#
def get_modified_id(id):
    # take care of characters that will mess up newick and/or fasttree
    row_id_disp = re.sub(r"\s", "_", id)
    row_id_disp = re.sub(r"\/", "%" + "/".encode("utf-8").hex(), row_id_disp)
    row_id_disp = re.sub(r"\\", "%" + "\\".encode("utf-8").hex(), row_id_disp)
    row_id_disp = re.sub(r"\(", "%" + "(".encode("utf-8").hex(), row_id_disp)
    row_id_disp = re.sub(r"\)", "%" + ")".encode("utf-8").hex(), row_id_disp)
    row_id_disp = re.sub(r"\[", "%" + "[".encode("utf-8").hex(), row_id_disp)
    row_id_disp = re.sub(r"\]", "%" + "]".encode("utf-8").hex(), row_id_disp)
    #row_id_disp = re.sub("\:", "%" + ":".encode("utf-8").hex(), row_id_disp)
    row_id_disp = re.sub(r"\:", "_", row_id_disp)
    row_id_disp = re.sub(r"\;", "_", row_id_disp)
    row_id_disp = re.sub(r"\|", "%" + ";".encode("utf-8").hex(), row_id_disp)
    row_id_disp = re.sub(r"\'", "", row_id_disp)

    return row_id_disp

# modify_ids_in_msa_for_newick ()
#
def modify_ids_in_msa_for_newick (gblocks_msa_input_path, data_dir, faa_file):
    gblocks_msa_modified_path = get_modified_msa_output_path (data_dir, faa_file)
    id_map = dict()
    outbuf = []

    with open (gblocks_msa_input_path,'r') as g_h:
        for g_line in g_h:
            g_line = g_line.strip()
            if g_line.startswith('>'):
                this_id = g_line.replace('>','')
                
                row_id_disp = get_modified_id(this_id)

                id_map[this_id] = row_id_disp
                
                new_line = g_line.replace(this_id,row_id_disp)
                outbuf.append(new_line)
            else:
                outbuf.append(g_line)
                

    with open (gblocks_msa_modified_path,'w') as g_h:
        g_h.write("\n".join(outbuf)+"\n")
        
    return (gblocks_msa_modified_path, id_map)

# run_muscle ()
#
def run_muscle (fam_seq_nr_input_path, data_dir, faa_file):

    print ("RUNNING MUSCLE...", flush=True)
    
    ### Construct the command
    #
    #  e.g. muscle -in <fasta_in> -out <fasta_out>
    #
    muscle_cmd = [MUSCLE_bin]

    # set the paths
    muscle_output_path = get_muscle_output_path (data_dir, faa_file)
    if os.path.exists(muscle_output_path) and \
       os.path.isfile(muscle_output_path) and \
       os.path.getsize(muscle_output_path) > 0:
        print("MUSCLE ALREADY RUN")
        return
    
    # build cmd
    # muscle_cmd.append('-align')
    muscle_cmd.append('-super5')
    muscle_cmd.append(fam_seq_nr_input_path)
    muscle_cmd.append('-output')
    muscle_cmd.append(muscle_output_path)

    # print(" ".join(muscle_cmd))

    # Run MUSCLE, capture output as it happens
    #
    p = subprocess.Popen(muscle_cmd, \
                         #cwd = data_dir, \
                         stdout = subprocess.PIPE, \
                         stderr = subprocess.STDOUT, \
                         shell = False)

    while True:
        line = p.stdout.readline().decode()
        if not line: break
        #self.log(console, line.replace('\n', ''))
        #print(line.replace('\n', ''))

    p.stdout.close()
    p.wait()
    #self.log(console, 'return code: ' + str(p.returncode))
    if p.returncode != 0:
        raise ValueError('Error running MUSCLE, return code: '+str(p.returncode) + '\n\n')

    return


# run_gblocks ()
#
def run_gblocks (muscle_path, data_dir, faa_file):

    print ("RUNNING GBLOCKS...", flush=True)

    ### Construct the command
    #
    #  e.g. muscle -in <fasta_in> -out <fasta_out> -maxiters <n> -maxhours <h>
    #
    gblocks_cmd = [GBLOCKS_bin]

    # make short id msa
    (short_id_msa_path, short_id_map) = make_short_id_msa (muscle_path)
    
    # set the paths
    short_id_msa_gblocks_out_path = short_id_msa_path+'-gb' # Gblocks makes it
    gblocks_output_path = get_gblocks_output_path (data_dir, faa_file)
    if os.path.exists(gblocks_output_path) and \
       os.path.isfile(gblocks_output_path) and \
       os.path.getsize(gblocks_output_path) > 0:
        print("GBLOCKS ALREADY RUN")
        return

    # run
    env = os.environ.copy()
    p = subprocess.Popen(gblocks_cmd, \
                         #cwd = self.scratch, \
                         stdin = subprocess.PIPE, \
                         stdout = subprocess.PIPE, \
                         stderr = subprocess.PIPE, \
                         shell = True, \
                         env = env)

    # write commands to process
    #
    #  for "0.5" gaps: cat "o\n<MSA_file>\nb\n5\ng\nm\nq\n" | Gblocks
    #  for "all" gaps: cat "o\n<MSA_file>\nb\n5\n5\ng\nm\nq\n" | Gblocks
    
    p.stdin.write(("o"+"\n").encode())  # open MSA file
    p.stdin.write((short_id_msa_path+"\n").encode())

    # set trim level
    default_trim_level = '1'
    p.stdin.write(("b"+"\n").encode())
    p.stdin.write(("5"+"\n").encode())  # set to "half" (default_trim_level 1)
    p.stdin.write(("m"+"\n").encode())

    # not doing flank
    default_flank = 0

    # not doing min_seqs_for_conserved
    default_min_seqs_for_conserved = 0

    # max_pos_nonconserved
    default_max_pos_nonconserved = 8
    p.stdin.write(("b"+"\n").encode())
    p.stdin.write(("3"+"\n").encode())
    p.stdin.write((str(default_max_pos_nonconserved)+"\n").encode())
    p.stdin.write(("m"+"\n").encode())

    # min block len
    default_min_block_len = 10
    p.stdin.write(("b"+"\n").encode())
    p.stdin.write(("4"+"\n").encode())
    p.stdin.write((str(default_min_block_len)+"\n").encode())
    p.stdin.write(("m"+"\n").encode())

    # run
    p.stdin.write(("g"+"\n").encode())  # get blocks
    p.stdin.write(("q"+"\n").encode())  # quit
    p.stdin.close()
    p.wait()

    while True:
        line = p.stdout.readline()
        #line = p.stderr.readline()
        if not line: break
        #print(line.replace('\n', ''))

    p.stdout.close()
    #p.stderr.close()
    p.wait()
    #self.log(console, 'return code: ' + str(p.returncode))
    #if p.returncode != 0:
    if p.returncode != 1:
        raise ValueError('Error running GBLOCKS, return code: '+str(p.returncode) + '\n\n')


    # make long id file without gaps and a single line fasta
    fix_gblocks_output (short_id_msa_gblocks_out_path,
                        gblocks_output_path,
                        short_id_map)
    
    # clean up intermediate tmp files
    os.remove(short_id_msa_path)
    os.remove(short_id_msa_gblocks_out_path)
    os.remove(short_id_msa_gblocks_out_path+'.htm')

    return


# run_fasttree ()
#
def run_fasttree (gblocks_input_path, data_dir, faa_file):

    print ("RUNNING FASTTREE...", flush=True)

    ### Construct the command
    #
    #  e.g. fasttree 
    #
    fasttree_cmd = [FASTTREE_bin]

    # set the paths
    tree_newick_output_path = get_tree_output_path (data_dir, faa_file)

    # fix ids to not clash with newick format
    (gblocks_msa_modified_path, id_map) = modify_ids_in_msa_for_newick (gblocks_input_path, data_dir, faa_file)
    id_map_path = get_id_map_path (data_dir, faa_file)

    if os.path.exists(tree_newick_output_path) and \
       os.path.isfile(tree_newick_output_path) and \
       os.path.getsize(tree_newick_output_path) > 0 and \
       os.path.exists(gblocks_msa_modified_path) and \
       os.path.isfile(gblocks_msa_modified_path) and \
       os.path.getsize(gblocks_msa_modified_path) > 0 and \
       os.path.exists(id_map_path) and \
       os.path.isfile(id_map_path) and \
       os.path.getsize(id_map_path) > 0:
        print("FASTTREE ALREADY RUN")
        return

    with open (id_map_path,'w') as id_h:
        for this_id in sorted(id_map.keys()):
            id_h.write("\t".join([this_id,id_map[this_id]])+"\n")

    # set defaults
    default_cat = 20
    
    fasttree_cmd.append("-nopr")
    #fasttree_cmd.append("-fastest")
    #fasttree_cmd.append("-pseudo")
    #fasttree_cmd.append("-intree")
    #fasttree_cmd.append(str(intree_newick_file_path))
    #fasttree_cmd.append("-noml")
    #fasttree_cmd.append("-nome")
    #fasttree_cmd.append("-nocat")
    fasttree_cmd.append("-cat")
    fasttree_cmd.append(str(default_cat))
    fasttree_cmd.append("-gamma")


    # This doesn't work for some reason
    #        fasttree_cmd.append('-out')
    #        fasttree_cmd.append(output_newick_file_path)

    # better (meaning it works) to write MSA to STDIN (below)
    fasttree_cmd.append(">")
    fasttree_cmd.append(str(tree_newick_output_path))
    
    
    # Run FASTREE, capture output as it happens
    #
    env = os.environ.copy()
    joined_fasttree_cmd = " ".join(fasttree_cmd)
    p = subprocess.Popen([joined_fasttree_cmd], \
                         #cwd = data_dir, \
                         stdin  = subprocess.PIPE,
                         stdout = subprocess.PIPE, \
                         stderr = subprocess.PIPE, \
                         shell = True,
                         env = env)

    with open (gblocks_msa_modified_path,'r') as g_h:
        for g_line in g_h:
            p.stdin.write(g_line.encode())
    p.stdin.close()
    p.wait()

    while True:
        line = p.stdout.readline().decode()
        if not line: break
        #self.log(console, line.replace('\n', ''))
        #print(line.replace('\n', ''))

    p.stdout.close()
    p.wait()
    #self.log(console, 'return code: ' + str(p.returncode))
    if p.returncode != 0:
        raise ValueError('Error running FASTTREE, return code: '+str(p.returncode) + '\n\n')

    # clean up
    os.remove(gblocks_msa_modified_path)

    return

# get_label_colors(bitscore_map, clade_names)
#
def get_label_colors(bitscore_map, clade_names):
    label_colors = dict()
    sm = None

    if bitscore_map: 
        cmap = plt.get_cmap('viridis')

        new_bitscore_map = {name: bitscore_map[name] for name in clade_names if name in bitscore_map.keys()}

        if new_bitscore_map:
            # normalize values
            max_score = max(new_bitscore_map.values())
            min_score = min(new_bitscore_map.values())
            if max_score - min_score != 0:
                norm_bitscore_map = {name: (score - min_score) / (max_score - min_score) for name, score in new_bitscore_map.items()}
            else:
                norm_bitscore_map = {name: 0 for name in new_bitscore_map.keys()}

            label_colors = {name: mcolors.to_hex(cmap(norm_score)) for name, norm_score in norm_bitscore_map.items()}
            
            sm = plt.cm.ScalarMappable(norm=Normalize(vmin=min_score, vmax=max_score), cmap=cmap)
    
    return label_colors, sm

# read_id_map(id_map_path)
#
def read_id_map(id_map_path):
    og_id_map = dict()
    with open(id_map_path) as file:
        lines = file.readlines()
        for line in lines:
            line = line.rstrip("\n")
            values = line.split("\t")
            og_id_map[values[1]] = values[0]
    
    return og_id_map

# vis_tree(fasttree_input_path, data_dir, faa_file, bitscore_map)
#
def vis_tree(fasttree_input_path, data_dir, faa_file, bitscore_map):
    print ("VISUALIZING...", flush=True)

    vis_output_path = get_phylo_vis_output_path(data_dir, faa_file)
    if os.path.exists(vis_output_path) and \
       os.path.isfile(vis_output_path) and \
       os.path.getsize(vis_output_path) > 0:
        print("VISUALIZING ALREADY RUN")
        return

    tree = Phylo.read(fasttree_input_path, "newick")

    clade_names = set()
    og_id_map = read_id_map(get_id_map_path (data_dir, faa_file))
    for clade in tree.get_terminals():
        clade.name = og_id_map[clade.name]
        clade_names.add(clade.name)

    depth = max(list(tree.depths().values()))
    width_factor = 5
    width = (depth * width_factor) + 60
    num_leaves = len(tree.get_terminals())
    height_per_leaf = 0.25
    height = max(10, num_leaves * height_per_leaf)
    
    label_colors, sm = get_label_colors(bitscore_map, clade_names)

    fig = plt.figure(figsize=(width, height), dpi=75)
    axes = fig.add_subplot(1, 1, 1)
    axes.set_facecolor("#a0a0a0")

    if sm:
        cbar = plt.colorbar(sm, ax=axes)
        cbar.set_label("Bitscore")

    Phylo.draw(tree, axes=axes, do_show=False, label_func=lambda x: x.name, label_colors=label_colors)
    plt.savefig(vis_output_path)
    plt.close()

# remove_repeat_seqs(out_dir, faa_file)
#
def remove_repeat_seqs(out_dir, faa_file):
    print ("DEPUPLICATING SEQUENCES...", flush=True)

    # set the paths
    output_path = get_fam_seq_nr_path (out_dir, faa_file)
    if os.path.exists(output_path) and \
       os.path.isfile(output_path) and \
       os.path.getsize(output_path) > 0:
        print("DEPUPLICATING SEQUENCES ALREADY RUN")
        return
    
    faa_dict = read_fasta_to_dict(faa_file)

    seen = set()
    remove = set()
    for name, seq in faa_dict.items():
        if seq in seen:
            remove.add(name)
        else:
            seen.add(seq)
    
    for name in remove:
        faa_dict.pop(name)

    write_dict_to_fasta(output_path, faa_dict)

# get_ref_positions(ref_seqs, faa_dict)
#
def get_ref_positions(ref_seqs, faa_dict):
    positions = defaultdict(set)
    for name in ref_seqs:
        seq = faa_dict[name]
        for pos, aa in enumerate(seq):
            if aa != "-":
                positions[name].add(pos)
    
    return positions

# filter_msa_by_coverage(muscle_output_path, out_dir, faa_file, ref_seq_map, filter_coverage)
#
def filter_msa_by_coverage(muscle_output_path, out_dir, faa_file, ref_seq_map, filter_coverage):
    print ("FILTERING MSA...", flush=True)
    
    # set the paths
    filtered_output_path = get_muscle_filtered_output_path (out_dir, faa_file)
    if os.path.exists(filtered_output_path) and \
       os.path.isfile(filtered_output_path) and \
       os.path.getsize(filtered_output_path) > 0:
        print("FILTERING MSA ALREADY RUN")
        return
    
    faa_dict = read_fasta_to_dict(muscle_output_path)

    ref_seqs = ref_seq_map[os.path.basename(faa_file).replace('.faa', '')]
    all_ref_positions = get_ref_positions(ref_seqs, faa_dict)

    keep = set(ref_seqs)
    for name, seq in faa_dict.items():
        if name not in ref_seqs:
            for ref in ref_seqs:
                ref_positions = all_ref_positions[ref]

                coverage_count = 0
                for pos in ref_positions:
                    if seq[pos] != "-":
                        coverage_count += 1

                if coverage_count >= len(ref_positions) * filter_coverage:
                    keep.add(name)
                    break

    faa_dict = {name: seq for name, seq in faa_dict.items() if name in keep}

    write_dict_to_fasta(filtered_output_path, faa_dict)

# get_filtered_out_seqs(filtered_faa_file, out_dir, faa_file)
#
def get_filtered_out_seqs(filtered_faa_file, out_dir, faa_file):
    filtered_list_output_path = get_filtered_list_output_path (out_dir, faa_file)
    if os.path.exists(filtered_list_output_path) and \
       os.path.isfile(filtered_list_output_path) and \
       os.path.getsize(filtered_list_output_path) > 0:
        return

    og_faa = read_fasta_to_dict(faa_file)
    filtered_faa = read_fasta_to_dict(filtered_faa_file)

    filtered_out_seqs = (set(og_faa.keys()) - set(filtered_faa.keys()))

    with open(filtered_list_output_path, "w") as file:
        for seq in filtered_out_seqs:
            file.write(f"{seq}\n")

# run_process(out_dir, faa_file, filter_coverage, refref_seqs_map_seq_map, bitscore_map)
#
def run_process(out_dir, faa_file, filter_coverage, ref_seqs_map, bitscore_map):
    if not more_than_one_seq_in(faa_file):
        print("SKIPPED - one or less sequence(s)")
        return
    
    ''' Note: not removing repeats because it may be the case that an ICTV ref seq is the same as a novel seq
    # remove repeat seqs
    remove_repeat_seqs(out_dir, faa_file)

    fam_seq_nr_input_path = get_fam_seq_nr_path (out_dir, faa_file)
    if not more_than_one_seq_in(fam_seq_nr_input_path):
        print("SKIPPED - one or less sequence(s)")
        return
    '''

    # run muscle
    run_muscle(faa_file, out_dir, faa_file)

    # filter to overlapping seqs
    muscle_output_path = get_muscle_output_path (out_dir, faa_file)
    filter_msa_by_coverage(muscle_output_path, out_dir, faa_file, ref_seqs_map, filter_coverage)

    muscle_filtered_output_path = get_muscle_filtered_output_path (out_dir, faa_file)
    if not more_than_one_seq_in(muscle_filtered_output_path):
        print("SKIPPED - one or less sequence(s)")
        return

    # run Gblocks
    run_gblocks(muscle_filtered_output_path, out_dir, faa_file)

    # run FastTree
    gblocks_output_path = get_gblocks_output_path (out_dir, faa_file)
    if not has_sequences(gblocks_output_path):
        print("SKIPPED - no sequences")
        return
    
    get_filtered_out_seqs(gblocks_output_path, out_dir, faa_file)
    
    run_fasttree(gblocks_output_path, out_dir, faa_file)

    fasttree_output_path = get_tree_output_path (out_dir, faa_file)
    vis_tree(fasttree_output_path, out_dir, faa_file, bitscore_map)        

    """
    # final file
    newick_tree_path = get_tree_output_path (args.datadir, args.faafile)
    print ("created {}".format(newick_tree_path))
    """

    pass

# has_sequences(muscle_output_path)
#
def has_sequences(muscle_output_path):
    with open(muscle_output_path, "r") as file:
        lines = file.readlines()
        for line in lines:
            line = line.replace("\n", "")
            if line and not line.startswith(">"):
                return True
    return False
    
# make_dir(path)
#
def make_dir(path):
    if not os.path.exists(path):
        os.mkdir(path)

# more_than_one_seq_in(faa_file)
#
def more_than_one_seq_in(faa_file):
    count = 0
    with open(faa_file, "r") as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith(">"):
                count += 1
            if count > 1:
                return True
    return False

# read_bitscore_map(bitscore_map_filepath)
#
def read_bitscore_map(bitscore_map_filepath):
    if bitscore_map_filepath:
        bitscore_map = dict()
        with open(bitscore_map_filepath) as file:
            lines = file.readlines()
            for line in lines:
                line = line.rstrip("\n")
                values = line.split("\t")
                bitscore_map[values[0]] = float(values[1])
        return bitscore_map
    return None

# combine_seqs(data_dir, ref_seq_dir, outdir)
#
def combine_seqs(data_dir, ref_seq_dir, outdir):
    out_faa_map = dict()
    ref_seqs_map = defaultdict(set)
    fam_file_map = defaultdict(list)

    for root, dirs, files in os.walk(data_dir):
        for file in files:
            faa_name = file.replace("-sliced.faa", "")
            fam_file_map[faa_name].append(os.path.join(root, file))

    for root, dirs, files in os.walk(ref_seq_dir):
        for file in files:
            faa_name = file.replace("-sliced.faa", "")
            if faa_name in fam_file_map.keys(): # check if faa is in previous set
                faa_filepath = os.path.join(root, file)
                fam_file_map[faa_name].append(faa_filepath)

                faa_dict = read_fasta_to_dict(faa_filepath)
                ref_seqs_map[faa_name] = set(faa_dict.keys())
    
    for faa, files in fam_file_map.items():
        faa_dir = os.path.join(outdir, faa)
        make_dir(faa_dir)
        out_faa_filepath = os.path.join(faa_dir, f"{faa}.faa")
        out_faa_map[faa] = out_faa_filepath

        with open(out_faa_filepath, "w") as w_file:
            for file in files:
                with open(file, "r") as r_file:
                    lines = r_file.readlines()
                    for line in lines:
                        w_file.write(line)
    
    return out_faa_map, ref_seqs_map

# main()
#
def main() -> int:
    args = getargs()

    bitscore_map = read_bitscore_map(args.bitscore_map)

    file_map, ref_seqs_map = combine_seqs(args.datadir, args.ref_seq_dir, args.outdir)

    for faa, fasta_file in file_map.items():
        print(f"\nProcessing {faa}...")
        fasta_outdir = os.path.join(args.outdir, faa)
        run_process(fasta_outdir, fasta_file, args.coverage, ref_seqs_map, bitscore_map)
    
    print ("DONE")
    return 0

# read_fasta_to_dict(filepath)
#
def read_fasta_to_dict(filepath):
    seq_map = dict()
    with open(filepath, "r") as file:
        lines = file.readlines()
        name = ""
        seq = ""
        for line in lines:
            line = line.strip("*\n")
            if line.startswith(">"):
                if name and seq:
                    seq_map[name] = seq
                name = line[1:]
                seq = ""
            else:
                seq += line
                
        if name and seq:
            seq_map[name] = seq
    
    return seq_map

# write_dict_to_fasta(filepath, seq_map)
#
def write_dict_to_fasta(filepath, seq_map):
    with open(filepath, "w") as file:
        for name, seq in seq_map.items():
            file.write(f">{name}\n")
            file.write(f"{seq}\n")

# exec()
#
if __name__ == '__main__':
    sys.exit(main())

