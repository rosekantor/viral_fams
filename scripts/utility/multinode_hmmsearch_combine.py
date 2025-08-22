#!/usr/bin/python3

import sys
import os
import argparse
import glob
import pandas as pd


# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="Combines all finished hmmsearch jobs")

    parser.add_argument("-e", "--execute-dir", help="Output directory frommultinode_hmmsearch_execute.py")
    parser.add_argument("-o", "--outdir", help="Output directory")
    
    args = parser.parse_args()
    args_pass = True

    if args.execute_dir is None:
        print ("must specify --{}\n".format('execute-dir'))
        args_pass = False
    elif not os.path.exists(args.execute_dir) or \
        not os.path.isdir(args.execute_dir):
        print ("--{} {} is not a dir\n".format('execute-dir', args.execute_dir))

    if args.outdir is None:
        print ("must specify --{}\n".format('outdir'))
        args_pass = False
    elif not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    elif not os.path.isdir(args.outdir):
        print ("--{} {} is not a dir\n".format('outdir', args.outdir))
        
    if not args_pass:
        parser.print_help()
        sys.exit (-1)

    return args


def check_files_complete(execute_dir):
    domtbl_files = sorted(glob.glob(os.path.join(execute_dir, "*_domtbl_out.txt")))

    if not domtbl_files:
        print("*_domtbl_out.txt files not found.")
        print("Exiting...")
        sys.exit(-1)

    incomplete = list()
    for filepath in domtbl_files:
        with open(filepath, "r") as file:
            last_line = file.readlines()[-1].strip()
            if last_line != "# [ok]":
                incomplete.append(filepath)
    
    if incomplete:
        print("The following outputs are incomplete:")
        for filepath in incomplete:
            print(filepath)
        print("Exiting...")
        sys.exit(-1)

    return domtbl_files


def combine_domtbl_files(outdir, domtbl_files):
    header = ""
    get_header = True
    footer = ""
    get_footer = True

    line_data = list()

    for filepath in domtbl_files:
        with open(filepath, "r") as file:
            for line in file:
                if line.startswith("#"):
                    if get_header:
                        header += line
                    elif get_footer:
                        footer += line
                else:
                    if get_header:
                        get_header = False

                    fields = line.split()
                    query = fields[3]
                    score = float(fields[7])

                    line_data.append((query, score, line))
        
        if get_footer:
            get_footer = False
    
    df = pd.DataFrame(line_data, columns=["query name", "score", "fullstring"]).sort_values(by=["query name", "score"], ascending=[True, False]).reset_index(drop=True)

    basename = os.path.basename("_".join(domtbl_files[0].split("_")[:-3]))
    outfile = os.path.join(outdir, f'{basename}_combined_domtbl_out.txt')
    with open(outfile, "w") as file:
        file.write(header)
        for fullstring in df["fullstring"]:
            file.write(fullstring)
        file.write(footer)


def main():
    args = getargs()

    domtbl_files = check_files_complete(args.execute_dir)

    combine_domtbl_files(args.outdir, domtbl_files)

    return 0


if __name__ == "__main__":
    main()