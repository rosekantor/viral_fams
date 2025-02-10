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

    parser.add_argument("-i", "--ingenomeidsfile", help="input genome ids list file")
    parser.add_argument("-g", "--genomedir", help="genome directory")
    parser.add_argument("-o", "--outtsvfile", help="output gene tsv table")
    
    args = parser.parse_args()
    args_pass = True

    if len(sys.argv) < 7:
        parser.print_help()
        sys.exit (-1)

    if args.ingenomeidsfile is None:
        print ("must specify --{}\n".format('inpdbchainids'))
        args_pass = False
    elif not os.path.exists(args.ingenomeidsfile) or \
         not os.path.isfile(args.ingenomeidsfile) or \
         not os.path.getsize(args.ingenomeidsfile) > 0:
        print ("--{} {} must exist and not be empty\n".format('ingenomeidsfile', args.ingenomeidsfile))
        args_pass = False

    if args.genomedir is None:
        print ("must specify --{}\n".format('genomedir'))
        args_pass = False
    elif not os.path.exists(args.genomedir) or \
         not os.path.isdir(args.genomedir):
        print ("--{} {} must exist and be a dir\n".format('genomedir', args.genomedir))
        args_pass = False
        
    if args.outtsvfile is None:
        print ("must specify --{}\n".format('outtsvfile'))
        args_pass = False

    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


# read_genome_ids()
#
def read_genome_ids (genome_ids_file):
    genome_ids = []

    with open (genome_ids_file, 'r') as gid_h:
        for gid_line in gid_h:
            gid_line = gid_line.strip()
            if not gid_line or gid_line.startswith('#'):
                continue
            gid = gid_line.split()[0]
            genome_ids.append(gid)

    return genome_ids


# find_surface_genes ()
#
def find_surface_genes (genomedir, genome_ids, surface_terms):
    surface_genes = dict()

    for gid in genome_ids:
        print ("READING GENOME {}".format(gid))
        protein_fasta_path = get_protein_fasta_path (genomedir, gid)
        protein_fasta_buf = get_protein_fasta_buf (protein_fasta_path)

        hits = dict()
        polyprotein = dict()
        for faa_line in protein_fasta_buf:

            if faa_line.startswith(">"):
                header_rec = faa_line.split()
                gene_id = header_rec[0].lstrip('>')
                gene_desc = ' '.join(header_rec[1:])
                gene_desc = re.sub(' \[[^\[]+\]$','',gene_desc)

                for term in surface_terms:
                    if term.lower() in gene_desc.lower():
                        hits[gene_id] = gene_desc
                        break
                if 'polyprotein' in gene_desc.lower():
                    polyprotein[gene_id] = gene_desc 
                    
        if hits:
            surface_genes[gid] = hits
        else:
            if polyprotein:
                surface_genes[gid] = polyprotein
            else:
                surface_genes[gid] = {'-':'-'}
            
    return surface_genes


# get_protein_fasta_path ()
#
def get_protein_fasta_path (genomedir, genome_id):
    protein_fasta_file = genome_id+'-protein.faa'
    return os.path.join(genomedir, genome_id, protein_fasta_file)
    

# get_protein_fasta_buf ()
#
def get_protein_fasta_buf (file_path):
    buf = []

    with open (file_path, 'r') as f_h:
        for f_line in f_h:
            f_line = f_line.rstrip("\n")
            buf.append(f_line)
            
    return buf
    

# write_merged_tsv_outfile ()
#
def write_merged_tsv_outfile (outtsvfile, genome_ids, surface_genes):
    outbuf = []

    # header
    row = ['genome_id', 'gene_ids', 'descs']
    outbuf.append("\t".join(row))

    for gid in genome_ids:
        row = [gid]
        gene_ids = []
        gene_descs = []
        for gene_id in surface_genes[gid].keys():
            gene_ids.append(gene_id)
            gene_descs.append(surface_genes[gid][gene_id])
        row.append(';'.join(gene_ids))
        row.append(';'.join(gene_descs))
        outbuf.append("\t".join(row))
                
    # write outfile
    with open(outtsvfile, 'w') as out_h:
        out_h.write("\n".join(outbuf)+"\n")
    
    return outtsvfile


# main()
#
def main() -> int:
    args = getargs()

    # config
    surface_terms = ['spike',
                     'glycoprotein',
                     'attachment',
                     'fusion',
                     'e1',
                     'e2'
    ]
    
    # get genomeids
    genome_ids = read_genome_ids (args.ingenomeidsfile)

    # read querys
    surface_genes = find_surface_genes (args.genomedir,
                                        genome_ids,
                                        surface_terms)

    # write merged tsv hit file
    out_file = write_merged_tsv_outfile (args.outtsvfile,
                                         genome_ids,
                                         surface_genes)
    
    print ("DONE")
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())

