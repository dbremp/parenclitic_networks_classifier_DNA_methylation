"""
This script includes the functions for the analysis of the topological features extracted from the parenclitic graphs.
These functions are used in topo_analysis.py and Classification.py
"""

import pandas as pd
from scipy.stats import ttest_ind
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, roc_auc_score
from sklearn.preprocessing import LabelEncoder


def data_formating(aggr, nonaggr):
    # this function brings the data in the right format for feature selection
    # from two separate dfs for aggr and non aggr it creates one including both
    # it also creates the y vector for the selector

    data = pd.concat([aggr, nonaggr], axis=0)
    aggr["aggressive"] = "yes"
    nonaggr["aggressive"] = "no"
    y = pd.concat([aggr["aggressive"], nonaggr["aggressive"]], axis=0)

    return data, y


def data_selection(data, y, k):
    # this function finds the k best features for data and y and returns the dataset reduced to these features
    # input: data, descriptor y, number of features to be selected k
    # output: data_new is a dataframe including the data of only the selected features,
    #         data_melt is dataframe with the selected data in "melted" format appropriate for boxplots
    #         features is a dataframe containing the selected features

    selector = SelectKBest(score_func=chi2, k=k)
    data_new = pd.DataFrame(selector.fit_transform(data, y))
    features_names = data.columns.values[selector.get_support()]

    data_new = data.loc[:, features_names]
    data_new["aggressive"] = y
    data_new["barcode"] = data_new.index
    data_melt = pd.melt(data_new, id_vars=["barcode", "aggressive"])

    return data_new, data_melt, features_names


def t_test(aggr, nonaggr, features_names):
    # this function executes a t-test between aggresive and non-aggressive on a list of features
    # input: aggr, nonaggr are dataframes including the values for aggressive and non-aggressive samples respectively
    #        features_names is a list of loci of interest
    # output: prints and returns the results of the t-test

    # t-test with unequal variance (Welch test)
    t = ttest_ind(aggr.loc[:, features_names], nonaggr.loc[:, features_names], axis=0, equal_var=False)
    # readable format
    t = pd.DataFrame(t, index=["statistic", "p-value"], columns=features_names)
    print(t)

    return t


def significant_loci_all_data(topofeat_aggr_tcga, topofeat_aggr_local, topofeat_nonaggr_tcga, topofeat_nonaggr_local,
                              features_names):
    # this function identifies the loci with significant differences in the topological feature values
    # taking into acocunt all samples (local and tcga - excluding the ones used for modelling)
    # input: topofeat_aggr_tcga, topofeat_aggr_local, topofeat_nonaggr_tcga, topofeat_nonaggr_local contain the values
    #        of the topological feature of interest for local/tcga aggressive/non-aggressive samples
    #        feature_names includes the loci of interest

    topofeat_aggr_all = pd.concat([topofeat_aggr_tcga, topofeat_aggr_local])
    topofeat_nonaggr_all = pd.concat([topofeat_nonaggr_tcga, topofeat_nonaggr_local])

    t_topofeat_all = t_test(topofeat_aggr_all, topofeat_nonaggr_all, features_names)
    sign_loci_topofeat_all = t_topofeat_all.columns[t_topofeat_all.loc["p-value"] < 0.05]

    sign_data_aggr = topofeat_aggr_all[sign_loci_topofeat_all]

    return sign_loci_topofeat_all, sign_data_aggr

def eval(y_true, y_pred):
    # this function evaluates the performance of a classifier and prints the results
    # input:
    # y_true an array containing the true labels
    # y_pred an array containing the labels predicted by the classifier

    print('\nAccuracy: %.2f'
          % accuracy_score(y_true, y_pred))
    print('Balanced accuracy: %.2f'
          % balanced_accuracy_score(y_true, y_pred))
    f1 = f1_score(y_true, y_pred, average=None)
    print("F1-score: no: %.2f, yes: %.2f" %(f1[0], f1[1]))
    lb = LabelEncoder()
    print("ROC-AUC: %.2f" % roc_auc_score(lb.fit_transform(y_true), lb.fit_transform(y_pred)))

