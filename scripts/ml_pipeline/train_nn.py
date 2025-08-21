#!/usr/bin/python3

import sys
import os
import argparse
import glob
import pandas as pd
from collections import defaultdict
from sklearn.neural_network import MLPClassifier
import numpy as np
import pickle

ORDERED_TAXA_RANKS = ["Domain", "Realm", "Kingdom", "Phylum", "Subphylum", "Class", "Order", "Suborder", "Family", "Subfamily", "Genus", "Species"]


class MultiOutputMLPClassifier:
    def __init__(self, n_outputs, classes):
        self._classifiers = [MLPClassifier(
                random_state=1, 
                max_iter=300
            ) for _ in range(n_outputs)]

        self._initialized = [False] * n_outputs
        
        self.classes = classes
        
    def partial_fit(self, X, y):
        for i, clf in enumerate(self._classifiers):
            if not self._initialized[i]:
                clf.partial_fit(X, y[:, i], classes=self.classes[i])
                self._initialized[i] = True
            else:
                clf.partial_fit(X, y[:, i])
    
    def predict(self, X):
        return np.column_stack([clf.predict(X) for clf in self._classifiers])

    def to_pickle(self, filepath):
        with open(filepath, 'wb') as file:
            pickle.dump(self, file)

def from_pickle(filepath):
    with open(filepath, 'rb') as file:
        return pickle.load(file)


# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-i", "--inputdata", help="")
    parser.add_argument("-e", "--encodings", help="")
    parser.add_argument("-c", "--chunksize", type=int, default=10000, help="")
    parser.add_argument("-o", "--outdir", help="")
    
    args = parser.parse_args()
    args_pass = True

    if args.inputdata is None:
        print ("must specify --{}\n".format('inputdata'))
        args_pass = False
    elif not os.path.exists(args.inputdata) or \
         not os.path.isfile(args.inputdata) or \
         not os.path.getsize(args.inputdata) > 0:
        print ("--{} {} is not a dir\n".format('inputdata', args.inputdata))
    
    if args.encodings is None:
        print ("must specify --{}\n".format('encodings'))
        args_pass = False
    elif not os.path.exists(args.encodings) or \
         not os.path.isfile(args.encodings) or \
         not os.path.getsize(args.encodings) > 0:
        print ("--{} {} is not a dir\n".format('encodings', args.encodings))

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


def read_rank_decodings(filepath):
    all_rank_decodings = defaultdict(dict)

    with open(filepath, "r") as file:
        lines = file.readlines()

    for line in lines[1:]:
        rank, lineage, encoding = line.strip().split("\t")
        all_rank_decodings[rank][int(encoding)] = lineage
    
    return all_rank_decodings


def main():
    args = getargs()

    all_rank_decodings = read_rank_decodings(args.encodings)

    all_columns = pd.read_csv(args.inputdata, sep="\t",  nrows=0)
    x_columns = [col for col in all_columns if col not in ["gene"] + ORDERED_TAXA_RANKS]
    y_columns = ORDERED_TAXA_RANKS
    classes = [np.array(list(all_rank_decodings[rank].keys())) for rank in ORDERED_TAXA_RANKS]

    clf = MultiOutputMLPClassifier(n_outputs = len(ORDERED_TAXA_RANKS), classes = classes)

    for i, chunk_df in enumerate(pd.read_csv(args.inputdata, sep="\t", chunksize=args.chunksize)):
        X = chunk_df[x_columns].to_numpy()
        y = chunk_df[y_columns].to_numpy()

        print(f"Training on chunk {i+1}...")
        clf.partial_fit(X, y)

    clf.to_pickle(os.path.join(args.outdir, "model.pkl"))

    print("Done")
    return 0


if __name__ == "__main__":
    main()