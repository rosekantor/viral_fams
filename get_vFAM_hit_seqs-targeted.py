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
    parser = argparse.ArgumentParser(description="get sequences for fam hits")

    parser.add_argument("-i", "--inputhmmhitsfile", help="merged fam hmmer hits tsv table file")
    parser.add_argument("-f", "--famids", help="fam ids (comma separated)")
    parser.add_argument("-s", "--sequencefile", help="source sequence file")
    parser.add_argument("-o", "--outlabel", help="output file")
    parser.add_argument("-c", "--coveragemin", help="hmm coverage min (def: 0.5)")
    
    args = parser.parse_args()
    args_pass = True

    if len(sys.argv) < 9:
        parser.print_help()
        sys.exit (-1)

    if args.famids is None:
        print ("must specify --{}\n".format('famids'))
        args_pass = False

    if args.inputhmmhitsfile is None:
        print ("must specify --{}\n".format('inputhmmhitsfile'))
        args_pass = False
    elif not os.path.exists(args.inputhmmhitsfile) or \
         not os.path.isfile(args.inputhmmhitsfile) or \
         not os.path.getsize(args.inputhmmhitsfile) > 0:
        print ("--{} {} must exist and not be empty\n".format('inputhmmhitsfile', args.inputhmmhitsfile))
        args_pass = False

    if args.sequencefile is None:
        print ("must specify --{}\n".format('sequencefile'))
        args_pass = False
    elif not os.path.exists(args.sequencefile) or \
         not os.path.isfile(args.sequencefile) or \
         not os.path.getsize(args.sequencefile) > 0:
        print ("--{} {} must exist and not be empty\n".format('sequencefile', args.sequencefile))
        args_pass = False

    if args.outlabel is None:
        print ("must specify --{}\n".format('outlabel'))
        args_pass = False

    if args.coveragemin is None:
        args.coveragemin = 0.5
    else:
        args.coveragemin = float(args.coveragemin)

    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


# get_genome_seq_path ()
#
def get_genome_seq_path (sequence_dir, genome_id):
    protein_seq_file = genome_id+'-protein.faa'
    genome_seq_path = os.path.join (sequence_dir, genome_id, protein_seq_file)
    return genome_seq_path


# get_fam_seq_path ()
#
def get_fam_seq_path (out_label, fam_id):
    dirpath = os.path.dirname(out_label)
    basename = os.path.basename(out_label)
    fam_seq_path = os.path.join(dirpath, fam_id+'-'+out_label+'-sliced.faa')
    return fam_seq_path


# read_genome_ids()
#
def read_genome_ids (genome_id_listfile):
    genome_ids = dict()
    
    with open (genome_id_listfile, 'r') as gid_h:
        for gid_line in gid_h:
            gid_line = gid_line.strip()
            if not gid_line or gid_line.startswith('#'):
                continue
            rec = gid_line.split("\t")
            genome_id = rec[0]
            genome_ids[genome_id] = True
                
    return genome_ids


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


# get_gene_seqs ()
#
def get_gene_seqs (genome_seq_file):
    headers = dict()
    gene_seqs = dict()
    last_gene_id = None
    seq = ''
    
    with open (genome_seq_file, 'r') as f_h:
        for faa_line in f_h:
            faa_line = faa_line.strip()
            if faa_line.startswith('>'):
                gene_id = faa_line.split()[0].replace('>','')
                headers[gene_id] = faa_line
                if seq:
                    gene_seqs[last_gene_id] = seq
                seq = ''
                last_gene_id = gene_id
            else:
                seq += faa_line
        if seq:
            gene_seqs[last_gene_id] = seq
            seq = ''
            last_gene_id = None

    return (gene_seqs, headers)


# read_genome_info()
#
def read_genome_info (genome_info_file):
    genome_info = dict()
    fields = []
    
    with open (genome_info_file, 'r') as ginfo_h:
        for ginfo_line in ginfo_h:
            ginfo_line = ginfo_line.strip()
            if not ginfo_line or ginfo_line.startswith('#'):
                continue
            if ginfo_line.startswith('rs_genome_id'):
                fields = ginfo_line.split("\t")
                continue
            rec = ginfo_line.split("\t")
            if not fields:
                raise ValueError ("need fields header in file {}".format(genome_info_file))
            genome_id = rec[0]
            genome_info[genome_id] = dict()
            for fi,field in enumerate(fields):
                genome_info[genome_id][field] = rec[fi]
                
    return (genome_info, fields)


# read_fam_hit_range ()
#
def read_fam_hit_range (fam_hits_file, fam_ids, coverage_min):
    fam_hits = dict()
    fam_hit_range = dict()
    fields = []
    #fi_lookup = dict()

    with open (fam_hits_file, 'r') as hmmer_h:
        for hmmer_line in hmmer_h:
            hmmer_line = hmmer_line.strip()
            if not hmmer_line or hmmer_line.startswith('#'):
                continue
            #if hmmer_line.startswith('genome_id'):
            if hmmer_line.startswith('target_name'):
                fields = hmmer_line.split("\t")
                #for fi,field in enumerate(fields):
                #    fi_lookup[field] = fi
                continue
            rec = hmmer_line.split("\t")
            if not fields:
                raise ValueError ("need fields header in file {}".format(fam_hits_file))
            hit = dict()
            for fi,field in enumerate(fields):
                hit[field] = rec[fi]
                
            fam_id = hit['target_name']
            #genome_id = hit['genome_id']
            gene_id = hit['query_name']

            #if fam_id not in fam_ids or genome_id not in genome_ids:
            if fam_id not in fam_ids:
                continue

#            if not valid_hmm_hit (coverage_min, hit['tlen'], hit['hmm_from'], hit['hmm_to']):
#                print ("{} to gene {} is too short.  coverage={}".format(fam_id, gene_id, ((float(hit['hmm_to'])-float(hit['hmm_from'])+1)/float(hit['tlen']))), flush=True)
#                continue
            
            if fam_id not in fam_hits:
                fam_hits[fam_id] = dict()
            #if genome_id not in fam_hits[fam_id]:
            #    fam_hits[fam_id][genome_id] = dict()
            #if gene_id not in fam_hits[fam_id][genome_id]:
            #    fam_hits[fam_id][genome_id][gene_id] = []
            #fam_hits[fam_id][genome_id][gene_id].append(hit)
            if gene_id not in fam_hits[fam_id]:
                fam_hits[fam_id][gene_id] = []
            fam_hits[fam_id][gene_id].append(hit)

    # reduce multi hits to range
    for fam_id in fam_hits.keys():
        fam_hit_range[fam_id] = dict()
        for gene_id in fam_hits[fam_id].keys():
            hmm_tlen = int(fam_hits[fam_id][gene_id][0]['tlen'])
            hmm_from = int(fam_hits[fam_id][gene_id][0]['hmm_from'])
            env_from = int(fam_hits[fam_id][gene_id][0]['env_from'])
            hmm_to = int(fam_hits[fam_id][gene_id][0]['hmm_to'])
            env_to = int(fam_hits[fam_id][gene_id][0]['env_to'])
            if len (fam_hits[fam_id][gene_id]) > 1:
                for hit in fam_hits[fam_id][gene_id][1:]:
                    if int(hit['hmm_from']) < hmm_from:
                        hmm_from = int(hit['hmm_from'])
                    if int(hit['env_from']) < env_from:
                        env_from = int(hit['env_from'])
                    if int(hit['hmm_to']) > hmm_to:
                        hmm_to = int(hit['hmm_to'])
                    if int(hit['env_to']) > env_to:
                        env_to = int(hit['env_to'])
                
            # validate coverage
            if not valid_hmm_hit (coverage_min, hmm_tlen, hmm_from, hmm_to):
                print ("{} to gene {} is too short.  coverage={}".format(fam_id, gene_id, ((float(hmm_to)-float(hmm_from)+1)/float(hmm_tlen))), flush=True)
                continue

            fam_hit_range[fam_id][gene_id] = dict()
            fam_hit_range[fam_id][gene_id]['env_from'] = env_from
            fam_hit_range[fam_id][gene_id]['env_to'] = env_to
            
    return fam_hit_range


# valid_hmm_hit()
#
def valid_hmm_hit (coverage_thresh, hmm_len, hmm_ali_from, hmm_ali_to):
    valid = False

    if float(int(hmm_ali_to)-int(hmm_ali_from)+1)/float(hmm_len) >= coverage_thresh:
        valid = True

    return valid

    
# write_fam_hit_sequences ()
#
def write_fam_hit_sequences (out_label,
                             sequence_file,
                             fam_hit_range):
    # build output
    for fam_id in sorted(fam_hit_range.keys()):
        print ("DOING {}".format(fam_id))
        outbuf = []
        

        # read gene sequences
        #genome_seq_path = get_genome_seq_path (sequence_dir, genome_id)
        genome_seq_path = sequence_file
        (gene_seqs, headers) = get_gene_seqs (genome_seq_path)

        # build sequence buf
        for gene_id in sorted(fam_hit_range[fam_id].keys()):
            print ("writing fasta for {}".format(gene_id),flush=True)  # DEBUG

            hit = fam_hit_range[fam_id][gene_id]
            #for hit in fam_hit_range[fam_id][gene_id]:
            env_from = hit['env_from']
            env_to = hit['env_to']
            seq_slice = gene_seqs[gene_id][int(env_from)-1:int(env_to)]
            #new_id = genome_id+'-f:'+gene_id+'-slice:'+env_from+'-'+env_to
            new_id = gene_id+'-slice:'+str(env_from)+'-'+str(env_to)
            mod_header = headers[gene_id].replace(gene_id,new_id)
            outbuf.append(mod_header)
            outbuf.append(seq_slice)
                
        # write outlabel
        fam_seq_path = get_fam_seq_path (out_label, fam_id)
        with open(fam_seq_path, 'w') as out_h:
            out_h.write("\n".join(outbuf)+"\n")
    
    return


# main()
#
def main() -> int:
    args = getargs()

    # get fam ids
    fam_ids = args.famids.split(',')

    # read fam hits to genomes and genes
    fam_hit_range = read_fam_hit_range (args.inputhmmhitsfile,
                                        fam_ids,
                                        args.coveragemin)
    
    # write sequences for each fam
    write_fam_hit_sequences (args.outlabel,
                             args.sequencefile,
                             fam_hit_range)
    
    print ("DONE")
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())

