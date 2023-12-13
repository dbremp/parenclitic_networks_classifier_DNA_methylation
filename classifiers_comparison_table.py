"""
This script applies a number of different classifiers to different topological features and
for each topological feature summarises the performance of all classifiers in a csv file.
"""

import os
import pandas as pd

from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import BaggingClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.linear_model import LogisticRegression

from classifiers_evaluation_functions import classifiers_evaluation

# point to the right directory
os.chdir("~/Data/Methylation")

# set loci number (n) and parenclitic networks threshold (t) for which you run the analysis to load the relevant data
n = 1000
t = 2
# base data used for the modelling of non-aggressive behaviour
base = "based_on_tcga_data"

# load data - csv file for each topological feature of interest
node_d = pd.read_csv("results/topo_features/" + base + "/tf_node_degrees_all_" + str(n) + "_" + str(t) + ".csv",
                     index_col=0)
btw_c = pd.read_csv("results/topo_features/" + base + "/tf_betweenness_centrality_all_" + str(n) + "_" + str(t) + ".csv",
                    index_col=0)
dg_c = pd.read_csv("results/topo_features/" + base + "/tf_degree_centrality_all_" + str(n) + "_" + str(t) + ".csv",
                   index_col=0)
ecc = pd.read_csv("results/topo_features/" + base + "/tf_eccentricity_all_" + str(n) + "_" + str(t) + ".csv",
                  index_col=0)
ecc = ecc.dropna()
eig_c = pd.read_csv("results/topo_features/" + base + "/tf_eigenvector_centrality_all_" + str(n) + "_" + str(t) + ".csv",
                    index_col=0)
sec_ord_c = pd.read_csv("results/topo_features/" + base + "/tf_second_order_centrality_all_" + str(n) + "_" + str(t) + ".csv",
                        index_col=0)
sec_ord_c = sec_ord_c.dropna()

# create y for local and tcga samples
tcga_samples = []
y_tcga = []
local_samples = []
y_local = []
for s in node_d.columns:
    if "TCGA" in s:
        tcga_samples = tcga_samples + [s]
        if "non" in s:
            y_tcga = y_tcga + ["no"]
        else:
            y_tcga = y_tcga + ["yes"]
    else:
        local_samples = local_samples + [s]
        if "non" in s:
            y_local = y_local + ["no"]
        else:
            y_local = y_local + ["yes"]


# test different features and classifiers to identify the best one
data = [node_d, btw_c, dg_c, ecc, eig_c, sec_ord_c]
data_names = ["node degree", "betweenness centrality", "degree centrality", "eccentricity", "eigenvector centrality", "second order centrality"]
cl_names = [
    "Logistic Regression",
    "Nearest Neighbors",
    "Decision Tree",
    "Neural Net",
    "AdaBoost",
    "Bagging KNN",
    "Bagging Logistic Regression"
]

random_state = 0
classifiers = [
    LogisticRegression(random_state=random_state),
    KNeighborsClassifier(3),
    DecisionTreeClassifier(max_depth=5),
    MLPClassifier(alpha=1, max_iter=1000),
    AdaBoostClassifier(),
    BaggingClassifier(KNeighborsClassifier(), max_samples=0.5, max_features=0.5, random_state=random_state),
    BaggingClassifier(LogisticRegression(random_state=random_state), max_samples=0.5, max_features=0.5)
]

for d, name in zip(data, data_names):
    classifiers_performance = classifiers_evaluation(d.loc[:, tcga_samples].transpose(),
                                                     d.loc[:, local_samples].transpose(), y_tcga, y_local,
                                                     cl_names, classifiers)
    classifiers_performance.to_csv("pnets_" + name + "_classifiers_comparison.csv")


print("done")
