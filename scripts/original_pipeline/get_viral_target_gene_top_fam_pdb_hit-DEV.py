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
    parser = argparse.ArgumentParser(description="combine tsvs for hmmer scans against viral gene families")

    parser.add_argument("-t", "--targetgenefile", help="viral target gene file")
    parser.add_argument("-f", "--famhitfile", help="viral fam hit file")
    parser.add_argument("-p", "--pdbblastdir", help="pdb blast dir")
    parser.add_argument("-o", "--outtsvfile", help="output family lineage and host reange tsv table")
    
    args = parser.parse_args()
    args_pass = True

    if len(sys.argv) < 9:
        parser.print_help()
        sys.exit (-1)

    if args.targetgenefile is None:
        print ("must specify --{}\n".format('targetgenefile'))
        args_pass = False
    elif not os.path.exists(args.targetgenefile) or \
         not os.path.isfile(args.targetgenefile) or \
         not os.path.getsize(args.targetgenefile) > 0:
        print ("--{} {} must exist and not be empty\n".format('targetgenefile', args.targetgenefile))
        args_pass = False

    if args.famhitfile is None:
        print ("must specify --{}\n".format('famhitfile'))
        args_pass = False
    elif not os.path.exists(args.famhitfile) or \
         not os.path.isfile(args.famhitfile) or \
         not os.path.getsize(args.famhitfile) > 0:
        print ("--{} {} must exist and not be empty\n".format('famhitfile', args.famhitfile))
        args_pass = False

    if args.pdbblastdir is None:
        print ("must specify --{}\n".format('pdbblastdir'))
        args_pass = False
    elif not os.path.exists(args.pdbblastdir) or \
         not os.path.isdir(args.pdbblastdir):
        print ("--{} {} must exist\n".format('pdbblastdir', args.pdbblastdir))
        args_pass = False

    if args.outtsvfile is None:
        print ("must specify --{}\n".format('outtsvfile'))
        args_pass = False

    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


# read_target_gene_file()
#
def read_target_gene_file (target_gene_file):
    target_genes = dict()
    fields = []
    
    with open (target_gene_file, 'r') as targ_h:
        for targ_line in targ_h:
            targ_line = targ_line.strip()
            if not targ_line or targ_line.startswith('#'):
                continue
            if targ_line.startswith("genome_id\t"):
                fields = targ_line.split("\t")
                continue
            rec = targ_line.split("\t")
            if not fields:
                raise ValueError ("need fields header in file {}".format(target_gene_file))
            this_targ = dict()
            for fi,field in enumerate(fields):
                this_targ[field] = rec[fi]
            genome_id = this_targ['genome_id']
            gene_ids = this_targ['gene_ids'].split(';')
            descs = this_targ['descs'].split(';')
            
            if genome_id not in target_genes:
                target_genes[genome_id] = dict()
            for gene_i,gene_id in enumerate(gene_ids): 
                target_genes[genome_id][gene_id] = descs[gene_i]
            
    return target_genes


# read_fam_top_hit()    REQUIRES HITS ARE ALREADY SORTED
#
def read_fam_top_hit (fam_hit_file, target_genes):
    fam_top_hit = dict()
    fields = []
    
    with open (fam_hit_file, 'r') as fhit_h:
        for fhit_line in fhit_h:
            fhit_line = fhit_line.strip()
            if not fhit_line or fhit_line.startswith('#'):
                continue
            if fhit_line.startswith("genome_id\t"):
                fields = fhit_line.split("\t")
                continue
            rec = fhit_line.split("\t")
            if not fields:
                raise ValueError ("need fields header in file {}".format(fam_hit_file))
            this_hit = dict()
            for fi,field in enumerate(fields):
                this_hit[field] = rec[fi]
            genome_id = this_hit['genome_id']
            gene_id = this_hit['query_name']
            fam_id = this_hit['target_name']

            if genome_id not in target_genes:
                continue
            if gene_id not in target_genes[genome_id]:
                continue
            
            if genome_id not in fam_top_hit:
                fam_top_hit[genome_id] = dict()
            if gene_id not in fam_top_hit[genome_id]:
                fam_top_hit[genome_id][gene_id] = fam_id
            
    return fam_top_hit


# read_pdb_top_hit()    REQUIRES HITS ARE ALREADY SORTED
#
def read_pdb_top_hit (pdb_blast_dir, target_genes):
    pdb_top_hit = dict()

    fields = 'query_acc subject_acc identity alnlen mismatches gap_opens q_beg q_end s_beg s_end e-value bitscore'.split()
    
    for genome_id in target_genes.keys():
        pdb_hit_path = get_pdb_hit_path (pdb_blast_dir, genome_id)
        pdb_hit_buf = get_pdb_hit_buf (pdb_hit_path)

        for pdb_hit_line in pdb_hit_buf:
            if not pdb_hit_line or pdb_hit_line.startswith('#'):
                continue
            rec = pdb_hit_line.split("\t")
            this_hit = dict() 
            for fi,field in enumerate(fields):
                this_hit[field] = rec[fi]
            gene_id = this_hit['query_acc']
            pdb_id = this_hit['subject_acc']
            pdb_id = fix_pdb_id (pdb_id)
            
            if gene_id not in target_genes[genome_id]:
                continue
            
            if genome_id not in pdb_top_hit:
                pdb_top_hit[genome_id] = dict()
            if gene_id not in pdb_top_hit[genome_id]:
                pdb_top_hit[genome_id][gene_id] = pdb_id
            
    return pdb_top_hit


# get_pdb_hit_path ()
#
def get_pdb_hit_path (pdb_dir, genome_id):
    pdb_hit_file = genome_id+'-protein-m7.pdb_blastp'
    return os.path.join(pdb_dir, pdb_hit_file)
    

# get_pdb_hit_buf ()
#
def get_pdb_hit_buf (file_path):
    buf = []

    with open (file_path, 'r') as f_h:
        for f_line in f_h:
            f_line = f_line.rstrip("\n")
            buf.append(f_line)
            
    return buf


# fix_pdb_id ()
#
def fix_pdb_id (pdb_id):
    new_pdb_id = None

    split_pdb_id = pdb_id.split('_')
    new_pdb_id = split_pdb_id[0]+'_'+','.join(list(split_pdb_id[1]))

    return new_pdb_id


# write_fam_hit_and_pdb_tsv_outfile ()
#
def write_fam_hit_and_pdb_tsv_outfile (outtsvfile,
                                       target_genes,
                                       fam_top_hit,
                                       pdb_top_hit):
    outbuf = []

    # header
    header_row = ['genome_id', 'gene_ids', 'top_hit_fams', 'top_hit_pdbs', 'descs']
    outbuf.append("\t".join(header_row))
    
    # genes
    for genome_id in sorted(target_genes.keys()):
        genes_buf = []
        descs_buf = []
        fam_top_hit_buf = []
        pdb_top_hit_buf = []
        
        for gene_id in sorted(target_genes[genome_id].keys()):
            genes_buf.append(gene_id)
            if gene_id == '-':
                descs_buf.append('-')
                fam_top_hit_buf.append('-')
                pdb_top_hit_buf.append('-')
            else:
                descs_buf.append(target_genes[genome_id][gene_id])
                if genome_id in fam_top_hit and gene_id in fam_top_hit[genome_id]:
                    fam_top_hit_buf.append(fam_top_hit[genome_id][gene_id])
                else:
                    fam_top_hit_buf.append('-')
                if genome_id in pdb_top_hit and gene_id in pdb_top_hit[genome_id]:
                    pdb_top_hit_buf.append(pdb_top_hit[genome_id][gene_id])
                else:
                    pdb_top_hit_buf.append('-')
                    
        outbuf.append("\t".join([genome_id,
                                 ';'.join(genes_buf),
                                 ';'.join(fam_top_hit_buf),
                                 ';'.join(pdb_top_hit_buf),
                                 ';'.join(descs_buf)
        ]))



    # write outfile
    if not outtsvfile:
        outtsvfile = 'STDOUT'
        print ("\n".join(outbuf))
    else:
        with open(outtsvfile, 'w') as out_h:
            out_h.write("\n".join(outbuf)+"\n")

    return outtsvfile


# main()
#
def main() -> int:
    args = getargs()

    # read genes
    target_genes = read_target_gene_file (args.targetgenefile)
    
    # read top hit to viral fam
    fam_top_hit = read_fam_top_hit (args.famhitfile, target_genes)

    # read top hit to pdb
    pdb_top_hit = read_pdb_top_hit (args.pdbblastdir, target_genes)
    
    # write lineage fams tsv outfile
    out_file = write_fam_hit_and_pdb_tsv_outfile (args.outtsvfile,
                                                  target_genes,
                                                  fam_top_hit,
                                                  pdb_top_hit)
    
    #print ("DONE")
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())

