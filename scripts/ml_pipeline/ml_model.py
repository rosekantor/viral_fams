# %%
from sklearn.multioutput import MultiOutputClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report
import pandas as pd
import numpy as np
import pickle
import os

# %%
ref_seq_table = "/p/lustre1/golez1/vog_231_ml_training_table_random_sample/training.tsv"

# %%
ORDERED_TAXA_RANKS = ["Domain", "Realm", "Kingdom", "Phylum", "Subphylum", "Class", "Order", "Suborder", "Family", "Subfamily", "Genus", "Species"]

# %%
header = pd.read_csv(ref_seq_table, sep="\t", nrows=0)
col_names = header.columns.tolist()

dtype_dict = {col_names[0]: str}
for col in col_names[1:-len(ORDERED_TAXA_RANKS)]:
    dtype_dict[col] = float
for col in col_names[-len(ORDERED_TAXA_RANKS):]:
    dtype_dict[col] = str

# %%
data_df = pd.read_csv(ref_seq_table, sep="\t", nrows=10000, dtype=dtype_dict)

x_columns = [col for col in data_df.columns if col not in ["gene"] + ORDERED_TAXA_RANKS]
y_columns = ORDERED_TAXA_RANKS[1:]  # remove Domain because it only has one class (Viruses)

X = data_df[x_columns].to_numpy()
y = data_df[y_columns].fillna("").to_numpy()

# %%
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state=1)

# %%
def test_model_accuracy(clf, X_test, y_test):
    y_pred = clf.predict(X_test)

    for i in range(len(y_columns)):
        test_col_arr = y_test[:, i]
        pred_col_arr = y_pred[:, i]
        
        print(f"Column: {y_columns[i]}")
        # print(f"Accuracy: {round(accuracy_score(test_col_arr, pred_col_arr), 2)}")
        print(classification_report(test_col_arr, pred_col_arr, zero_division=0))
        print()

# %%
def save_model(filepath, model):
    with open(filepath, 'wb') as file:
        pickle.dump(model, file)

def load_model(filepath):
    with open(filepath, 'rb') as file:
        return pickle.load(file)

# %% [markdown]
# Neural Network

# %%

from sklearn.neural_network import MLPClassifier

clf = MultiOutputClassifier(MLPClassifier(random_state=1, max_iter=300)).fit(X_train, y_train)

# %%
save_model(os.path.join("/p/lustre1/golez1/vog_231_ml_training_table_random_sample", "neural_network.pkl"), clf)

# %%
test_model_accuracy(clf, X_test, y_test)

# %% [markdown]
# Logistic Regression

# %%
from sklearn.linear_model import LogisticRegression

lr = MultiOutputClassifier(LogisticRegression(solver="saga", max_iter=1000)).fit(X_train, y_train)

# %%
save_model(os.path.join("/p/lustre1/golez1/vog_231_ml_training_table_random_sample", "logistic_regression.pkl"), lr)

# %%
test_model_accuracy(lr, X_test, y_test)

# %% [markdown]
# KNN

# %%
from sklearn.neighbors import KNeighborsClassifier

knn = MultiOutputClassifier(KNeighborsClassifier(n_neighbors=3)).fit(X_train, y_train)

# %%
save_model(os.path.join("/p/lustre1/golez1/vog_231_ml_training_table_random_sample", "knn.pkl"), knn)

# %%
test_model_accuracy(knn, X_test, y_test)

# %% [markdown]
# Random Forest

# %%
from sklearn.ensemble import RandomForestClassifier

forest = MultiOutputClassifier(RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42)).fit(X_train, y_train)

# %%
save_model(os.path.join("/p/lustre1/golez1/vog_231_ml_training_table_random_sample", "random_forest.pkl"), forest)

# %%
test_model_accuracy(forest, X_test, y_test)

# %% [markdown]
# Naive Bayes

# %%
from sklearn.naive_bayes import GaussianNB
nb = MultiOutputClassifier(GaussianNB()).fit(X_train, y_train)

# %%
save_model(os.path.join("/p/lustre1/golez1/vog_231_ml_training_table_random_sample", "naive_bayes.pkl"), nb)

# %%
test_model_accuracy(nb, X_test, y_test)

# %% [markdown]
# Gradient Boosting

# %%
from sklearn.ensemble import GradientBoostingClassifier

gb = MultiOutputClassifier(GradientBoostingClassifier(n_estimators=100, learning_rate=1.0, max_depth=1, random_state=42)).fit(X_train, y_train)

# %%
save_model(os.path.join("/p/lustre1/golez1/vog_231_ml_training_table_random_sample", "gradient_boosting.pkl"), gb)

# %%
test_model_accuracy(gb, X_test, y_test)


