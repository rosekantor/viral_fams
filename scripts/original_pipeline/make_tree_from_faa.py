#!/usr/bin/python3
'''
Copyright 2024 Dylan Chivian

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/gpl-3.0.txt>.
'''

import sys
import os
import argparse
import re
import subprocess
import shutil


# software run
#
#HOME_DIR = '~chivian1'
HOME_DIR = '/g/g16/chivian1'
FASTA_UTIL_DIR = os.path.join (HOME_DIR, 'src', 'fasta_util')
SOFTWARE_DIR = os.path.join (HOME_DIR, 'work', 'software')
MUSCLE_bin = os.path.join (SOFTWARE_DIR, 'MUSCLE_5.3', 'muscle-linux-x86.v5.3')
GBLOCKS_bin = os.path.join (SOFTWARE_DIR, 'Gblocks_0.91b', 'Gblocks')
FASTTREE_bin = os.path.join (SOFTWARE_DIR, 'FastTree_2.1.11', 'FastTreeDbl-linux-x86.v2.1.11')
REMOVE_REPEAT_SEQS_bin = os.path.join (FASTA_UTIL_DIR, 'remove_repeat_seqs.py')
FILTER_MSA_BY_REF_REGION_bin = os.path.join (FASTA_UTIL_DIR, 'filter_msa_by_ref_region_coverage.py')


# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="make tree from sequences for fam hits")

    parser.add_argument("-f", "--faafile", help="protein fasta sequence file")
    parser.add_argument("-d", "--datadir", help="data dir")
    parser.add_argument("-r", "--referenceid", help="reference id of gene with overlap region (optional)")
    parser.add_argument("-b", "--beginpos", help="begin pos of overlap region (optional)")
    parser.add_argument("-e", "--endpos", help="end pos of overlap region (optional)")
    parser.add_argument("-c", "--coverage", help="fraction of overlap region required covered [0.0,1.0] (def: 0.5)")

    
    args = parser.parse_args()
    args_pass = True

    if len(sys.argv) < 5:
        parser.print_help()
        sys.exit (-1)

    if args.faafile is None:
        print ("must specify --{}\n".format('faafile'))
        args_pass = False
    elif not os.path.exists(args.faafile) or \
         not os.path.isfile(args.faafile) or \
         not os.path.getsize(args.faafile) > 0:
        print ("--{} {} must exist and not be empty\n".format('faafile', args.faafile))
        args_pass = False

    if args.datadir is None:
        print ("must specify --{}\n".format('datadir'))
        args_pass = False
    elif not os.path.exists(args.datadir) or \
         not os.path.isdir(args.datadir):
        print ("--{} {} must exist\n".format('datadir', args.datadir))
        args_pass = False

    if not args.coverage:
        args.coverage = '0.5'
        
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
    muscle_filtered_path = os.path.join (data_dir, os.path.basename(faa_file).replace('.faa','-nr-muscle_filt.afa'))
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
    gblocks_path = os.path.join (data_dir, os.path.basename(faa_file).replace('.faa','-nr-muscle_filt-gblocks.afa'))
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
    mod_gblocks_path = os.path.join (data_dir, os.path.basename(faa_file).replace('.faa','-nr-muscle_filt-gblocks-modified.afa'))
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
    newick_tree_path = os.path.join (data_dir, os.path.basename(faa_file).replace('.faa','-nr-muscle_filt-gblocks-fasttree.newick'))
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
    id_map_path = os.path.join (data_dir, os.path.basename(faa_file).replace('.faa','-nr-muscle_filt-gblocks-fasttree.id_map'))
    return id_map_path


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
                this_id = g_line.split()[0].replace('>','')
                
                # take care of characters that will mess up newick and/or fasttree
                row_id_disp = re.sub("\s", "_", this_id)
                row_id_disp = re.sub("\/", "%" + "/".encode("utf-8").hex(), row_id_disp)
                row_id_disp = re.sub(r"\\", "%" + "\\".encode("utf-8").hex(), row_id_disp)
                row_id_disp = re.sub("\(", "%" + "(".encode("utf-8").hex(), row_id_disp)
                row_id_disp = re.sub("\)", "%" + ")".encode("utf-8").hex(), row_id_disp)
                row_id_disp = re.sub("\[", "%" + "[".encode("utf-8").hex(), row_id_disp)
                row_id_disp = re.sub("\]", "%" + "]".encode("utf-8").hex(), row_id_disp)
                #row_id_disp = re.sub("\:", "%" + ":".encode("utf-8").hex(), row_id_disp)
                row_id_disp = re.sub("\:", "_", row_id_disp)
                row_id_disp = re.sub("\;", "%" + ";".encode("utf-8").hex(), row_id_disp)
                row_id_disp = re.sub("\|", "%" + ";".encode("utf-8").hex(), row_id_disp)
                id_map[this_id] = row_id_disp
                
                new_line = g_line.replace(this_id,row_id_disp)
                outbuf.append(new_line)
            else:
                outbuf.append(g_line)
                

    with open (gblocks_msa_modified_path,'w') as g_h:
        g_h.write("\n".join(outbuf)+"\n")
        
    return (gblocks_msa_modified_path, id_map)


# run_uniq_seq ()
#
def run_uniq_seq (fam_seq_input_path, data_dir, faa_file): 
    print ("RUNNING UNIQ SEQ...", flush=True)

    rm_repeat_seqs_cmd = [REMOVE_REPEAT_SEQS_bin]

    # check for necessary files
    if not os.path.isfile(REMOVE_REPEAT_SEQS_bin):
        raise ValueError("no such file '"+MUSCLE_bin+"'")
    if not os.path.isfile(fam_seq_input_path):
        raise ValueError("no such file '"+fam_seq_input_path+"'")
    elif not os.path.getsize(fam_seq_input_path) > 0:
        raise ValueError("empty file '"+fam_seq_input_path+"'")

    # set the paths
    uniq_seq_path = get_fam_seq_nr_path (data_dir, faa_file)
    if os.path.exists(uniq_seq_path) and \
       os.path.isfile(uniq_seq_path) and \
       os.path.getsize(uniq_seq_path) > 0:
        print("UNIQ SEQ ALREADY RUN")
        return
    
    # build cmd
    rm_repeat_seqs_cmd.append('-i')
    rm_repeat_seqs_cmd.append(fam_seq_input_path)
    rm_repeat_seqs_cmd.append('-o')
    rm_repeat_seqs_cmd.append(uniq_seq_path)

    # Run, capture output as it happens
    #
    p = subprocess.Popen(rm_repeat_seqs_cmd, \
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
        raise ValueError('Error running REMOVE_REPEAT_SEQS, return code: '+str(p.returncode) + '\n\n')

    return


# run_filter_overlap ()
#
def run_filter_overlap (muscle_output_path, data_dir, faa_file,
                        reference_id, begin_pos, end_pos, min_coverage):
    print ("RUNNING FILTER MSA BY OVERLAP...", flush=True)

    filter_msa_cmd  = [FILTER_MSA_BY_REF_REGION_bin]

    # check for necessary files
    if not os.path.isfile(FILTER_MSA_BY_REF_REGION_bin):
        raise ValueError("no such file '"+MUSCLE_bin+"'")
    if not os.path.isfile(muscle_output_path):
        raise ValueError("no such file '"+muscle_output_path+"'")
    elif not os.path.getsize(muscle_output_path) > 0:
        raise ValueError("empty file '"+muscle_output_path+"'")

    # set the paths
    muscle_filtered_output_path = get_muscle_filtered_output_path (data_dir, faa_file)
    if os.path.exists(muscle_filtered_output_path) and \
       os.path.isfile(muscle_filtered_output_path) and \
       os.path.getsize(muscle_filtered_output_path) > 0:
        print("FILTER MSA ALREADY RUN")
        return

    if not reference_id:
        print ("NO MSA FILTER REQUESTED. JUST COPYING MSA TO FILT MSA")
        shutil.copy (muscle_output_path, muscle_filtered_output_path)
        return
    
    # build cmd
    filter_msa_cmd.append('-i')
    filter_msa_cmd.append(muscle_output_path)
    filter_msa_cmd.append('-o')
    filter_msa_cmd.append(muscle_filtered_output_path)
    filter_msa_cmd.append('-r')
    filter_msa_cmd.append(reference_id)
    filter_msa_cmd.append('-b')
    filter_msa_cmd.append(begin_pos)
    filter_msa_cmd.append('-e')
    filter_msa_cmd.append(end_pos)
    filter_msa_cmd.append('-c')
    filter_msa_cmd.append(min_coverage)
    

    # Run, capture output as it happens
    #
    p = subprocess.Popen(filter_msa_cmd, \
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
        raise ValueError('Error running FILTER_MSA, return code: '+str(p.returncode) + '\n\n')

    return


# run_muscle ()
#
def run_muscle (fam_seq_input_path, data_dir, faa_file):

    print ("RUNNING MUSCLE...", flush=True)
    
    ### Construct the command
    #
    #  e.g. muscle -in <fasta_in> -out <fasta_out> -maxiters <n> -maxhours <h>
    #
    muscle_cmd = [MUSCLE_bin]

    # check for necessary files
    if not os.path.isfile(MUSCLE_bin):
        raise ValueError("no such file '"+MUSCLE_bin+"'")
    if not os.path.isfile(fam_seq_input_path):
        raise ValueError("no such file '"+fam_seq_input_path+"'")
    elif not os.path.getsize(fam_seq_input_path) > 0:
        raise ValueError("empty file '"+fam_seq_input_path+"'")

    # set the paths
    muscle_output_path = get_muscle_output_path (data_dir, faa_file)
    if os.path.exists(muscle_output_path) and \
       os.path.isfile(muscle_output_path) and \
       os.path.getsize(muscle_output_path) > 0:
        print("MUSCLE ALREADY RUN")
        return
    
    # build cmd
    muscle_cmd.append('-align')
    muscle_cmd.append(fam_seq_input_path)
    muscle_cmd.append('-output')
    muscle_cmd.append(muscle_output_path)

    # options
    muscle_cmd.append('-maxiters')
    muscle_cmd.append(str(10))
    #muscle_cmd.append('-maxhours')
    #muscle_cmd.append(str(2))

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
def run_gblocks (muscle_msa_path, data_dir, faa_file):

    print ("RUNNING GBLOCKS...", flush=True)

    ### Construct the command
    #
    #  e.g. muscle -in <fasta_in> -out <fasta_out> -maxiters <n> -maxhours <h>
    #
    gblocks_cmd = [GBLOCKS_bin]

    # check for necessary files
    if not os.path.isfile(GBLOCKS_bin):
        raise ValueError("no such file '"+GBLOCKS_bin+"'")
    if not os.path.isfile(muscle_msa_path):
        raise ValueError("no such file '"+muscle_msa_path+"'")
    elif not os.path.getsize(muscle_msa_path) > 0:
        raise ValueError("empty file '"+muscle_msa_path+"'")

    # make short id msa
    (short_id_msa_path, short_id_map) = make_short_id_msa (muscle_msa_path)
    
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
def run_fasttree (gblocks_msa_input_path, data_dir, faa_file):

    print ("RUNNING FASTTREE...", flush=True)

    ### Construct the command
    #
    #  e.g. fasttree 
    #
    fasttree_cmd = [FASTTREE_bin]

    # check for necessary files
    if not os.path.isfile(FASTTREE_bin):
        raise ValueError("no such file '"+FASTREE_bin+"'")
    if not os.path.isfile(gblocks_msa_input_path):
        raise ValueError("no such file '"+gblocks_msa_input_path+"'")
    elif not os.path.getsize(gblocks_msa_input_path) > 0:
        raise ValueError("empty file '"+gblocks_msa_input_path+"'")

    # set the paths
    tree_newick_output_path = get_tree_output_path (data_dir, faa_file)

    # fix ids to not clash with newick format
    (gblocks_msa_modified_path, id_map) = modify_ids_in_msa_for_newick (gblocks_msa_input_path, data_dir, faa_file)
    id_map_path = get_id_map_path (data_dir, faa_file)
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


# main()
#
def main() -> int:
    args = getargs()

    # remove redundant seqs
    fam_seq_input_path = get_fam_seq_path (args.datadir, args.faafile)
    run_uniq_seq (fam_seq_input_path, args.datadir, args.faafile)

    # run muscle
    fam_seq_nr_input_path = get_fam_seq_nr_path (args.datadir, args.faafile)
    run_muscle (fam_seq_nr_input_path, args.datadir, args.faafile)

    # filter to overlapping seqs
    muscle_output_path = get_muscle_output_path (args.datadir, args.faafile)
    run_filter_overlap (muscle_output_path, args.datadir, args.faafile,
                        args.referenceid, args.beginpos, args.endpos, args.coverage)
    
    # run Gblocks
    muscle_filtered_output_path = get_muscle_filtered_output_path (args.datadir, args.faafile)
    run_gblocks (muscle_filtered_output_path, args.datadir, args.faafile)

    # run FastTree
    gblocks_output_path = get_gblocks_output_path (args.datadir, args.faafile)
    run_fasttree (gblocks_output_path, args.datadir, args.faafile)

    # final file
    newick_tree_path = get_tree_output_path (args.datadir, args.faafile)
    print ("created {}".format(newick_tree_path))
        
        
    print ("DONE")
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())

