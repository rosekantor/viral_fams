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


# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="create ITOL annotation files for trees")

    parser.add_argument("-m", "--mapfile", help="input id map file")
    parser.add_argument("-f", "--fastafile", help="input fasta file")
    parser.add_argument("-l", "--lineagefile", help="list of ids by lineage (e.g. human")
    parser.add_argument("-p", "--pdbhitfile", help="input pdb hit file")
    parser.add_argument("-o", "--outdir", help="output dir for itol annot files")
    
    args = parser.parse_args()
    args_pass = True

    if len(sys.argv) < 11:
        parser.print_help()
        sys.exit (-1)

    if args.mapfile is None:
        print ("must specify --{}\n".format('mapfile'))
        args_pass = False
    elif not os.path.exists(args.mapfile) or \
         not os.path.isfile(args.mapfile) or \
         not os.path.getsize(args.mapfile) > 0:
        print ("--{} {} must exist and not be empty\n".format('mapfile', args.mapfile))
        args_pass = False

    if args.fastafile is None:
        print ("must specify --{}\n".format('fastafile'))
        args_pass = False
    elif not os.path.exists(args.fastafile) or \
         not os.path.isfile(args.fastafile) or \
         not os.path.getsize(args.fastafile) > 0:
        print ("--{} {} must exist and not be empty\n".format('fastafile', args.fastafile))
        args_pass = False

    if args.lineagefile is None:
        print ("must specify --{}\n".format('lineagefile'))
        args_pass = False
    elif not os.path.exists(args.lineagefile) or \
         not os.path.isfile(args.lineagefile) or \
         not os.path.getsize(args.lineagefile) > 0:
        print ("--{} {} must exist and not be empty\n".format('lineagefile', args.lineagefile))
        args_pass = False

    if args.pdbhitfile is None:
        print ("must specify --{}\n".format('pdbhitfile'))
        args_pass = False
    elif not os.path.exists(args.pdbhitfile) or \
         not os.path.isfile(args.pdbhitfile) or \
         not os.path.getsize(args.pdbhitfile) > 0:
        print ("--{} {} must exist and not be empty\n".format('pdbhitfile', args.pdbhitfile))
        args_pass = False

    if args.outdir is None:
        print ("must specify --{}\n".format('outdir'))
        args_pass = False

    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


# read_id_map()
#
def read_id_map (map_file):
    id_map = dict()
    with open (map_file,'r') as m_h:
        for line in m_h:
            line = line.rstrip()
            (orig_id, leaf_id) = line.split()
            id_map[orig_id] = leaf_id
    return id_map


# read_species_names()
#
def read_species_names (fasta_file):
    sp_names = dict()
    with open (fasta_file,'r') as m_h:
        for line in m_h:
            line = line.rstrip()
            if line.startswith('>'):
                sp_name = re.sub('^.*\[','',line)
                sp_name = re.sub('\].*$','',sp_name)
                orig_id = line.split()[0].replace('>','')
                sp_names[orig_id] = sp_name
    return sp_names


# read_lineage_ids()
#
def read_lineage_ids (lineage_file):
    lineage_ids = dict()
    with open (lineage_file,'r') as m_h:
        for line in m_h:
            line = line.rstrip()
            genome_id = line.split()[0]
            lineage_ids[genome_id] = True
    return lineage_ids


# read_pdb_hits()
#
def read_pdb_hits (pdb_hits_file):
    pdb_hits = dict()
    fields = []
    with open (pdb_hits_file,'r') as m_h:
        for line in m_h:
            line = line.rstrip()
            if line.startswith ('v_genome_id'):
                fields = line.split("\t")
                continue
            row = line.split("\t")
            rec = dict()
            for fi,field in enumerate(fields):
                rec[field] = row[fi]
            combo_id = rec['v_genome_id']+'-f:'+rec['v_gene_id']
            if rec['pdb_hit_ident'] != '-':
                pdb_hits[combo_id] = rec['pdb_hit_ident']
    return pdb_hits


# write_id_map()
#
def write_id_map (outfile, id_map):
    outbuf = []
    for this_id in sorted(id_map.keys()):
        outbuf.append("\t".join([this_id,id_map[this_id]]))

    with open (outfile,'w') as o_h:
        o_h.write("\n".join(outbuf)+"\n")
        
    return


# write_leaf_labels ()
#
def write_leaf_labels (outdir, id_map, sp_names):
    outbuf = []
    outfile = os.path.join(outdir, 'annot_leaf_labels.itol')

    outbuf.append('LABELS')
    outbuf.append('SEPARATOR TAB')
    outbuf.append('DATA')
    
    for orig_id in sorted(id_map.keys()):
        leaf_id = id_map[orig_id]
        sp_name = sp_names[orig_id]
        if sp_name.startswith('gbkey'):
            sp_name = re.sub ('-slice:.*$','',orig_id)
        outbuf.append("\t".join([leaf_id,sp_name]))

    with open (outfile,'w') as o_h:
        o_h.write("\n".join(outbuf)+"\n")
    return

# write_label_colors ()
#
def write_label_colors (outdir, id_map, human_ids):
    outbuf = []
    outfile = os.path.join(outdir, 'annot_label_colors.itol')

    outbuf.append('TREE_COLORS')
    outbuf.append('SEPARATOR TAB')
    outbuf.append('DATA')

    for orig_id in sorted(id_map.keys()):
        leaf_id = id_map[orig_id]
        genome_id = re.sub('-f:.*$','', orig_id)
        if genome_id in human_ids:
            outbuf.append("\t".join([leaf_id,'label','#0000ff', 'bold', '1']))
        elif genome_id.startswith('O'):
            outbuf.append("\t".join([leaf_id,'label','#00cc00', 'bold', '1']))
            
    with open (outfile,'w') as o_h:
        o_h.write("\n".join(outbuf)+"\n")
    return


# write_pdb_barplot ()
#
def write_pdb_barplot (outdir, id_map, pdb_hits):
    outbuf = []
    outfile = os.path.join(outdir, 'annot_pdb_bars.itol')

    outbuf.append('DATASET_SIMPLEBAR')
    outbuf.append('SEPARATOR SPACE')
    outbuf.append('DATASET_LABEL PDB_HIT')
    outbuf.append('COLOR #000000')
    outbuf.append('DATA')

    for orig_id in sorted(id_map.keys()):
        base_id = re.sub('-slice:.*$','',orig_id)
        leaf_id = id_map[orig_id]
        if base_id in pdb_hits:
            #print ("COMBO_ID: {} {}".format(combo_id, rec['pdb_hit_ident']))
            outbuf.append(" ".join([leaf_id,str(pdb_hits[base_id])]))
            
    with open (outfile,'w') as o_h:
        o_h.write("\n".join(outbuf)+"\n")
    return

    

# main()
#
def main() -> int:
    args = getargs()

    id_map = read_id_map (args.mapfile)
    
    species_names = read_species_names (args.fastafile)

    human_ids = read_lineage_ids (args.lineagefile)
    
    pdb_hits = read_pdb_hits (args.pdbhitfile) 

    # write itol annot files
    write_leaf_labels (args.outdir, id_map, species_names)
    write_label_colors (args.outdir, id_map, human_ids)
    write_pdb_barplot (args.outdir, id_map, pdb_hits)
    
    #write_id_map (args.mappingfile, short_id_map)

    
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())

