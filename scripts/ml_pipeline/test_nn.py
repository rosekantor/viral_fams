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
from sklearn.metrics import accuracy_score, classification_report

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
            pickle.dump(self._classifiers, file)

def from_pickle(filepath):
    with open(filepath, 'rb') as file:
        return pickle.load(file)


# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-i", "--inputdata", help="")
    parser.add_argument("-m", "--model", help="")
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

    if args.model is None:
        print ("must specify --{}\n".format('model'))
        args_pass = False
    elif not os.path.exists(args.model) or \
         not os.path.isfile(args.model) or \
         not os.path.getsize(args.model) > 0:
        print ("--{} {} is not a dir\n".format('model', args.model))
    
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


def decode_predictions(y_pred, all_rank_decodings):
    max_len = 0
    for rank in ORDERED_TAXA_RANKS:
        for decoding in all_rank_decodings[rank].values():
            if len(decoding) > max_len:
                max_len = len(decoding)

    y_pred_decoded = np.zeros((y_pred.shape[0], y_pred.shape[1]), dtype=f'<U{max_len}')
    for i, rank in enumerate(ORDERED_TAXA_RANKS):
        for j, encoding in enumerate(y_pred[:, i]):
            y_pred_decoded[j, i] = all_rank_decodings[rank][encoding]
    
    return y_pred_decoded


def test_model_accuracy(clf_list, all_rank_decodings, genes, X_test, y_test):
    rank_accuracies = dict()
    predicted_data = list()

    # y_pred = clf.predict(X_test)
    y_pred = np.column_stack([clf.predict(X_test) for clf in clf_list])

    for i, rank in enumerate(ORDERED_TAXA_RANKS):
        test_col_arr = y_test[:, i]
        pred_col_arr = y_pred[:, i]

        # print(f"{round(len([name for name in test_col_arr if name not in  all_rank_decodings[rank].values()]) / len(test_col_arr) * 100, 2)} % not in possible classifications")
        accuracy = round(accuracy_score(test_col_arr, pred_col_arr), 2)

        # print(f"Column: {ORDERED_TAXA_RANKS[i]}")
        # print(f"Accuracy: {accuracy}")
        # print(classification_report(test_col_arr, pred_col_arr, zero_division=0))
        # print()

        rank_accuracies[rank] = accuracy
    
    y_pred_decoded = decode_predictions(y_pred, all_rank_decodings)
    
    for i in range(y_pred_decoded.shape[0]):
        info = [genes[i]]
        for j in range(y_pred_decoded.shape[1]):
            info.append(y_pred_decoded[i, j])
        
        predicted_data.append(tuple(info))
    
    # print(rank_accuracies)
    # print(predicted_data)
    
    return rank_accuracies, predicted_data


def write_dict_to_file(outfile, map):
    with open(outfile, "w") as file:
        for key, val in map.items():
            file.write(f"{key}\t{val}\n")


def main():
    args = getargs()

    all_rank_decodings = read_rank_decodings(args.encodings)

    all_columns = pd.read_csv(args.inputdata, sep="\t",  nrows=0)
    x_columns = [col for col in all_columns if col not in ["gene"] + ORDERED_TAXA_RANKS]
    y_columns = ORDERED_TAXA_RANKS
    # classes = [np.array(list(all_rank_decodings[rank].keys())) for rank in ORDERED_TAXA_RANKS]

    clf_list = from_pickle(args.model)

    rank_accuracies = {rank: 0 for rank in ORDERED_TAXA_RANKS}
    predictions = list()

    num_chunks = 0
    for i, chunk_df in enumerate(pd.read_csv(args.inputdata, sep="\t", chunksize=args.chunksize)):
        genes = chunk_df["gene"].to_numpy()
        X = chunk_df[x_columns].to_numpy()
        y = chunk_df[y_columns].to_numpy()

        print(f"Testing chunk {i+1}...")

        chunk_rank_accuracies, chunk_predictions = test_model_accuracy(clf_list, all_rank_decodings, genes, X, y)

        for rank in ORDERED_TAXA_RANKS:
            rank_accuracies[rank] += chunk_rank_accuracies[rank]

        predictions += chunk_predictions

        num_chunks += 1
    
    for rank in ORDERED_TAXA_RANKS:
        rank_accuracies[rank] = rank_accuracies[rank] / num_chunks

    write_dict_to_file(os.path.join(args.outdir, "accuracy.tsv"), rank_accuracies)
    
    pd.DataFrame(
        predictions, 
        columns=["gene"] + ORDERED_TAXA_RANKS
    ).to_csv(
        os.path.join(args.outdir, "predictions.tsv"),
        sep="\t",
        index=False
    )

    print("Done")
    return 0


if __name__ == "__main__":
    main()