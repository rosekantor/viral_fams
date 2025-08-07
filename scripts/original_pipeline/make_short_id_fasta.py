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


# software run
#
#HOME_DIR = '~chivian1'
#HOME_DIR = '/g/g16/chivian1'
#SOFTWARE_DIR = os.path.join (HOME_DIR, 'work', 'software')
#MUSCLE_bin = os.path.join (SOFTWARE_DIR, 'MUSCLE_5.3', 'muscle-linux-x86.v5.3')
#GBLOCKS_bin = os.path.join (SOFTWARE_DIR, 'Gblocks_0.91b', 'Gblocks')
#FASTTREE_bin = os.path.join (SOFTWARE_DIR, 'FastTree_2.1.11', 'FastTreeDbl-linux-x86.v2.1.11')


# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="get sequences for fam hits")

    parser.add_argument("-i", "--infastafile", help="input fasta file")
    parser.add_argument("-o", "--outfastafile", help="output fasta file with short ids")
    parser.add_argument("-m", "--mappingfile", help="output mapping file of short to long id")
    
    args = parser.parse_args()
    args_pass = True

    if len(sys.argv) < 5:
        parser.print_help()
        sys.exit (-1)

    if args.infastafile is None:
        print ("must specify --{}\n".format('infastafile'))
        args_pass = False
    elif not os.path.exists(args.infastafile) or \
         not os.path.isfile(args.infastafile) or \
         not os.path.getsize(args.infastafile) > 0:
        print ("--{} {} must exist and not be empty\n".format('infastafile', args.infastafile))
        args_pass = False

    if args.outfastafile is None:
        print ("must specify --{}\n".format('outfastafile'))
        args_pass = False

    if args.mappingfile is None:
        print ("must specify --{}\n".format('mappingfile'))
        args_pass = False


    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


# make_short_id_msa ()
#
def make_short_id_msa (muscle_msa_path, short_id_msa_path):
    short_id_map = dict()
    #short_id_msa_path = muscle_msa_path+'-short_ids'

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


# write_id_map()
#
def write_id_map (outfile, id_map):
    outbuf = []
    for this_id in sorted(id_map.keys()):
        outbuf.append("\t".join([this_id,id_map[this_id]]))

    with open (outfile,'w') as o_h:
        o_h.write("\n".join(outbuf)+"\n")
        
    return


# main()
#
def main() -> int:
    args = getargs()

    # make short id msa
    (short_id_msa_path, short_id_map) = make_short_id_msa (args.infastafile, args.outfastafile)

    write_id_map (args.mappingfile, short_id_map)

    
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())

