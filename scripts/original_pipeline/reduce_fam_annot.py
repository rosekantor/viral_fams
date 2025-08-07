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
    parser = argparse.ArgumentParser(description="reformat and reduce viral family annotation file")

    parser.add_argument("-i", "--infamannotfile", help="input family annot file")
    parser.add_argument("-o", "--outtsvfile", help="output reduced family annot file")
    
    args = parser.parse_args()
    args_pass = True

    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit (-1)

    if args.infamannotfile is None:
        print ("must specify --{}\n".format('inpdbchainids'))
        args_pass = False
    elif not os.path.exists(args.infamannotfile) or \
         not os.path.isfile(args.infamannotfile) or \
         not os.path.getsize(args.infamannotfile) > 0:
        print ("--{} {} must exist and not be empty\n".format('infamannotfile', args.infamannotfile))
        args_pass = False

    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


# read_fam_annots
#
def read_fam_annots (fam_annot_file):
    fam_annots = dict()

    with open (fam_annot_file, 'r') as fa_h:
        for fa_line in fa_h:
            fa_line = fa_line.strip()
            if not fa_line or fa_line.startswith('#'):
                continue
            [fam_rec,annot_rec] = fa_line.split("\t")
            if fam_rec == '-':
                continue
            if ';' not in fam_rec:
                if fam_rec not in fam_annots:
                    fam_annots[fam_rec] = dict()
                fam_annots[fam_rec][annot_rec] = True
            else:
                fams = fam_rec.split(';')
                annots = annot_rec.split(';')
                for rec_i,fam in enumerate(fams):
                    if fam == '-':
                        continue
                    if fam not in fam_annots:
                        fam_annots[fam] = dict()
                    fam_annots[fam][annots[rec_i]] = True
                
    return fam_annots


# write_tsv_outfile()
#
def write_tsv_outfile (outtsvfile, fam_annots):
    outbuf = []

    for fam in sorted (fam_annots.keys()):
        for annot in sorted(fam_annots[fam]):
            outbuf.append("\t".join([fam,annot]))

    # write outfile
    if outtsvfile:
        with open(outtsvfile, 'w') as out_h:
            out_h.write("\n".join(outbuf)+"\n")
    else:
        outtsvfile = 'STDOUT'
        print ("\n".join(outbuf))
        
    return outtsvfile


# main()
#
def main() -> int:
    args = getargs()

    fam_annots = read_fam_annots (args.infamannotfile)

    # write tsv out file
    out_file = write_tsv_outfile (args.outtsvfile, fam_annots)
    
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())

