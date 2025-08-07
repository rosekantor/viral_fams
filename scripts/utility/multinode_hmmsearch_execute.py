#!/usr/bin/python3

import sys
import os
import argparse
import subprocess


# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-f", "--faa", help="")
    parser.add_argument("-p", "--hmm", help="")
    parser.add_argument("-b", "--basename", help="")
    parser.add_argument("-n", "--nodes", type=int, default=1, help="")
    parser.add_argument("-c", "--cpus", type=int, default=100, help="")
    parser.add_argument("-o", "--outdir", help="")
    
    args = parser.parse_args()
    args_pass = True

    if args.faa is None:
        print ("must specify --{}\n".format('faa'))
        args_pass = False
    elif not os.path.exists(args.faa) or \
         not os.path.isfile(args.faa) or \
         not os.path.getsize(args.faa) > 0:
        print ("--{} {} must exist and not be empty\n".format('faa', args.faa))
        args_pass = False

    if args.hmm is None:
        print ("must specify --{}\n".format('hmm'))
        args_pass = False
    elif not os.path.exists(args.hmm) or \
         not os.path.isfile(args.hmm) or \
         not os.path.getsize(args.hmm) > 0:
        print ("--{} {} must exist and not be empty\n".format('hmm', args.hmm))
        args_pass = False
    
    if args.basename is None:
        print ("must specify --{}\n".format('basename'))
        args_pass = False
    
    if args.nodes is None:
        print ("must specify --{}\n".format('nodes'))
        args_pass = False

    if args.cpus is None:
        print ("must specify --{}\n".format('cpus'))
        args_pass = False

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


# read_fasta_to_dict()
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


# make_dir()
#
def make_dir(path):
    if not os.path.exists(path):
        os.mkdir(path)


# write_dict_to_fasta()
#
def write_dict_to_fasta(filepath, seq_map):
    with open(filepath, "w") as file:
        for name, seq in seq_map.items():
            file.write(f">{name}\n")
            file.write(f"{seq}\n")


# write_split_fasta()
#
def write_split_fasta(split_faa_outdir, seq_map, basename, num_nodes):
    split_faas = list()

    prot_ids = list(seq_map.keys())
    num_sequences = len(prot_ids)

    part_size = num_sequences // num_nodes
    bins_with_one_more = num_sequences % num_nodes

    stop = 0
    for i in range(num_nodes):
        start = stop
        if i < bins_with_one_more:
            stop += part_size + 1
        elif i < num_nodes - 1:
            stop += part_size
        else:
            stop = num_sequences

        partial_seq_map = {prot_id: seq_map[prot_id] for prot_id in prot_ids[start:stop]}

        basename_no_extension, extension = os.path.splitext(basename)
        partial_faa_out = os.path.join(split_faa_outdir, f"{basename_no_extension}_part_{i}{extension}")
        write_dict_to_fasta(partial_faa_out, partial_seq_map)

        split_faas.append(partial_faa_out)
    
    return split_faas


# write_bash_scripts()
#
def write_bash_scripts(outdir, logs_outdir, bash_script_outdir, split_faas, hmm_path, basename, cpus):
    bash_scripts = list()
    for i, faa_path in enumerate(split_faas):
        bash_script_out = os.path.join(bash_script_outdir, f"{basename}_job_{i}.sh")
        with open(bash_script_out, "w") as file:
            script = (
                f"#!/bin/bash\n"
                f"\n"
                f"#SBATCH -t 00-24:00:00\n"
                f"#SBATCH -N 1\n"
                f"#SBATCH -o {os.path.join(logs_outdir, f'hmmsearch_job_{i}.out')}\n"
                f"#SBATCH -e {os.path.join(logs_outdir, f'hmmsearch_job_{i}.err')}\n"
                f"#SBATCH --job-name {basename}_{i}\n"
                f"\n"
                f"srun hmmsearch \\\n"
                f"-o {os.path.join(outdir, f'{basename}_{i}_outfile.txt')} \\\n"
                f"--tblout {os.path.join(outdir, f'{basename}_{i}_tbl_out.txt')} \\\n"
                f"--domtblout {os.path.join(outdir, f'{basename}_{i}_domtbl_out.txt')} \\\n"
                f"--notextw \\\n"
                f"--cpu {cpus} \\\n"
                f"{hmm_path} {faa_path}"
            )

            file.write(script)

        bash_scripts.append(bash_script_out)

    return bash_scripts


# cancel_jobs()
#
def cancel_jobs(job_ids):
    for job_id in job_ids:
        try:
            result = subprocess.run(
                ['scancel', job_id],
                capture_output=True,
                text=True,
                check=True
            )
            print(f"Cancelled job {job_id}")

        except subprocess.CalledProcessError as e:
            print(f"Error cancelling job {job_id}: {e.stderr.strip()}")

# execute_scripts()
#
def execute_scripts(bash_scripts):
    submitted_job_ids = list()

    for script_path in bash_scripts:
        try:
            result = subprocess.run(
                ['sbatch', script_path],
                capture_output=True,
                text=True,
                check=True
            )
            print(f"Submitted {script_path}: {result.stdout.strip()}")

            job_id = str(result.stdout).split()[-1]
            submitted_job_ids.append(job_id)

        except subprocess.CalledProcessError as e:
            print(f"Error submitting {script_path}: {e.stderr.strip()}")
            print(f"Cancelling all other jobs...")
            cancel_jobs(submitted_job_ids)

            return list()

    return submitted_job_ids


def main():
    args = getargs()

    seq_map = read_fasta_to_dict(args.faa)

    split_faa_outdir = os.path.join(args.outdir, "split_faas")
    make_dir(split_faa_outdir)
    split_faas = write_split_fasta(split_faa_outdir, seq_map, os.path.basename(args.faa), args.nodes)

    bash_script_outdir = os.path.join(args.outdir, "bash_scripts")
    make_dir(bash_script_outdir)
    logs_outdir = os.path.join(args.outdir, "logs")
    make_dir(logs_outdir)
    bash_scripts = write_bash_scripts(args.outdir, logs_outdir, bash_script_outdir, split_faas, args.hmm, args.basename, args.cpus)

    submitted_job_ids = execute_scripts(bash_scripts)
     
    return 0


if __name__ == "__main__":
    main()