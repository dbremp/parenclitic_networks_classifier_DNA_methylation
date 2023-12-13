"""
This script includes the Logistic Regression classifier and the evaluation of its performance.
It also contains the identification of the most important loci for the classification,
the bar plot of their coefficient values and their annotation.
"""

import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn import metrics

from topo_analysis_functions import eval

# point to the right directory
os.chdir("~Data/Methylation")

# set the loci number (n) and threshold (t) for which you run the analysis in order to load the relevant data
n = 1000
t = 2
base = "based_on_tcga_data"

# load data
# betweenness centrality
node_d = pd.read_csv("results/topo_features/" + base + "/tf_node_degrees_all_" + str(n) + "_" + str(t) + ".csv", index_col=0)
anno = pd.read_csv("~/Data/Methylation/anno_minfi.csv", index_col=0)

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

# Logistic Regression classifier
# evaluate performance over 50 different random states
print("\nlogistic regression")
acc_scores = []
bacc_scores = []
roc_auc_scores = []
for random_state in range(50):
    classifier = LogisticRegression(random_state=random_state)
    classifier.fit(node_d.loc[:, tcga_samples].transpose(), pd.get_dummies(y_tcga)["yes"])
    y_pred = classifier.predict(node_d.loc[:, local_samples].transpose())
    acc_scores = acc_scores + [metrics.accuracy_score(pd.get_dummies(y_local)["yes"], y_pred)]
    bacc_scores = bacc_scores + [metrics.balanced_accuracy_score(pd.get_dummies(y_local)["yes"], y_pred)]
    roc_auc_scores = roc_auc_scores + [metrics.roc_auc_score(pd.get_dummies(y_local)["yes"], y_pred)]

print("accuracy:" + str(np.average(acc_scores)) + " +- " + str(np.std(acc_scores)))
print("balanced accuracy:" + str(np.average(bacc_scores)) + " +- " + str(np.std(bacc_scores)))
print("ROC-AUC:" + str(np.average(roc_auc_scores)) + " +- " + str(np.std(roc_auc_scores)))

# Logistic Regression classification
classifier = LogisticRegression(random_state=0)
classifier.fit(node_d.loc[:, tcga_samples].transpose(), pd.get_dummies(y_tcga)["yes"])
y_pred = classifier.predict(node_d.loc[:, local_samples].transpose())
eval(pd.get_dummies(y_local)["yes"], y_pred)

# identify the loci that play an important role for the classification
# sort coefficients
important_features = classifier.coef_
important_loci_df = pd.DataFrame(important_features.transpose(), index=node_d.index, columns=["coefficient"])
important_loci_df.sort_values(by="coefficient", ascending=False).to_csv("./results/classification/LogisticRegression/"
                                                                        "coefficients_node_d.csv")

# pick highest 5% abs(coef)
important_loci_df.quantile(q=0.975) # 0.001842
important_loci_df.quantile(q=0.025) # 0.001782
important_loci_df = important_loci_df[abs(important_loci_df["coefficient"]) > 0.0018]
sorted_important_loci_df = important_loci_df.sort_values(by="coefficient", ascending=False)
sorted_important_loci_df["cg_locus"] = sorted_important_loci_df.index

# create barplot showing the coefficient values of the important loci
plt.figure(figsize=(8,10))
ax = sns.barplot(x=sorted_important_loci_df.coefficient, y=sorted_important_loci_df.cg_locus,
                 hue=(sorted_important_loci_df.coefficient > 0))
ax.get_legend().remove()
# plt.show(block=True)
plt.savefig("./results/pnet_plots/coefficients_barchart.png")
plt.savefig("./results/pnet_plots/coefficients_barchart.pdf")

# include annotation of the important loci and save as csv table
anno_important_loci_df = anno.loc[sorted_important_loci_df.index]

sorted_important_loci_df = pd.merge(sorted_important_loci_df,
                                    anno_important_loci_df[["chr","Probe_rs","CpG_rs", "Relation_to_Island",
                                                            "UCSC_RefGene_Name","DMR", "X450k_Enhancer",
                                                            "Regulatory_Feature_Group","GencodeBasicV12_NAME",
                                                            "OpenChromatin_NAME", "TFBS_NAME"]],
                                    left_on="0", right_on="0")
sorted_important_loci_df["UCSC_RefGene_Name"] = [gene.split(";")[0]
                                                 for gene in
                                                 list(map(str, sorted_important_loci_df["UCSC_RefGene_Name"]))]
sorted_important_loci_df["GencodeBasicV12_NAME"] = [gene.split(";")[0]
                                                    for gene in
                                                    list(map(str, sorted_important_loci_df["GencodeBasicV12_NAME"]))]
sorted_important_loci_df.to_csv("./results/classification/LogisticRegression/important_loci_anno_node_d.csv")
