"""
This script:
* selects the k=1000 best features for the data and the phenotype label (y)
* find loci with significant differences in the topological features values between aggressive and non-aggressive
(considering both the TCGA and the Array Express samples, excluding the samples used for modelling)
* creates the boxplots for each topological features in the loci with significant differences
"""

import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from topo_analysis_functions import data_formating, data_selection, t_test, significant_loci_all_data

# point to the right directory
os.chdir("~/Data/Methylation/")
# set the loci number (n) and threshold (t) for which you run the analysis in order to load the relevant data
n = 1000
t = 2

base = "based_on_tcga_data"

# number of features I want to select
k = 1000

# load data
# betweenness centrality
btw_c_aggr = pd.read_csv("results/topo_features/" + base + "/betweenness_centrality_aggressive_" + str(n) + "_" +
                         str(t) + ".csv", index_col=0)
btw_c_nonaggr = pd.read_csv("results/topo_features/" + base + "/betweenness_centrality_non_aggressive_" + str(n) + "_" +
                            str(t) + ".csv", index_col=0)
btw_c_aggr_tcga = pd.read_csv("results/topo_features/" + base + "/betweenness_centrality_aggressive_tcga_" + str(n) +
                              "_" + str(t) + ".csv", index_col=0)
btw_c_nonaggr_tcga = pd.read_csv("results/topo_features/" + base + "/betweenness_centrality_non_aggressive_tcga_" +
                                 str(n) + "_" + str(t) + ".csv", index_col=0)

# degree centrality
dg_c_aggr = pd.read_csv("results/topo_features/" + base + "/degree_centrality_aggressive_" + str(n) + "_" + str(t) +
                        ".csv", index_col=0)
dg_c_nonaggr = pd.read_csv("results/topo_features/" + base + "/degree_centrality_non_aggressive_" + str(n) + "_" +
                           str(t) + ".csv", index_col=0)
dg_c_aggr_tcga = pd.read_csv("results/topo_features/" + base + "/degree_centrality_aggressive_tcga_" + str(n) + "_" +
                             str(t) + ".csv", index_col=0)
dg_c_nonaggr_tcga = pd.read_csv("results/topo_features/" + base + "/degree_centrality_non_aggressive_tcga_" + str(n) +
                                "_" + str(t) + ".csv", index_col=0)

# node degree
node_d_aggr = pd.read_csv("results/topo_features/" + base + "/node_degrees_aggressive_" + str(n) + "_" + str(t) +
                          ".csv", index_col=0)
node_d_nonaggr = pd.read_csv("results/topo_features/" + base + "/node_degrees_non_aggressive_" + str(n) + "_" +
                             str(t) + ".csv", index_col=0)
node_d_aggr_tcga = pd.read_csv("results/topo_features/" + base + "/node_degrees_aggressive_tcga_" + str(n) + "_" +
                               str(t) + ".csv", index_col=0)
node_d_nonaggr_tcga = pd.read_csv("results/topo_features/" + base + "/node_degrees_non_aggressive_tcga_" + str(n) +
                                  "_" + str(t) + ".csv", index_col=0)

# bring data to the right format
btw_c, y = data_formating(btw_c_aggr, btw_c_nonaggr)
btw_c_tcga, y_tcga = data_formating(btw_c_aggr_tcga, btw_c_nonaggr_tcga)
dg_c, y2 = data_formating(dg_c_aggr, dg_c_nonaggr)
dg_c_tcga, y2_tcga = data_formating(dg_c_aggr_tcga, dg_c_nonaggr_tcga)
node_d, y3 = data_formating(node_d_aggr, node_d_nonaggr)
node_d_tcga, y3_tcga = data_formating(node_d_aggr_tcga,node_d_nonaggr_tcga)

# sanity check
# it must be y == y2 == y2 and y_tcga == y2_tcga == y3_tcga
# this is important because only y and y_tcga are used for the rest of the analysis
print("sanity check:")
if ((y==y2).all() & (y2==y3).all() & (y_tcga==y2_tcga).all() & (y2_tcga==y3_tcga).all() & (y==y2).all() &
        (y2==y3).all() & (y_tcga==y2_tcga).all() & (y2_tcga==y3_tcga).all()):
    print("all good so far")
else:
    print("ATTENTION: The descriptors do not coincide! ")


# feature selection and formatting data for boxplots
btw_c_new, btw_c_melt, features_names_btw = data_selection(btw_c, y, k)
btw_c_new_tcga, btw_c_melt_tcga, features_names_btw_tcga = data_selection(btw_c_tcga, y_tcga, k)

dg_c_new, dg_c_melt, features_names_dg = data_selection(dg_c, y, k)
dg_c_new_tcga, dg_c_melt_tcga, features_names_dg_tcga = data_selection(dg_c_tcga, y_tcga, k)

node_d_new, node_d_melt, features_names_node_d = data_selection(node_d, y, k)
node_d_new_tcga, node_d_melt_tcga, features_names_node_d_tcga = data_selection(node_d_tcga, y_tcga, k)

# t-test
t_btw = t_test(btw_c_aggr, btw_c_nonaggr, features_names_btw)
t_btw_tcga = t_test(btw_c_aggr_tcga, btw_c_nonaggr_tcga, features_names_btw_tcga)
t_dg = t_test(dg_c_aggr, dg_c_nonaggr, features_names_dg)
t_dg_tcga = t_test(dg_c_aggr_tcga, dg_c_nonaggr_tcga, features_names_dg_tcga)
t_node_d = t_test(node_d_aggr, node_d_nonaggr, features_names_node_d)
t_node_d_tcga = t_test(node_d_aggr_tcga, node_d_nonaggr_tcga, features_names_node_d_tcga)

# find significantly different loci btw aggr and non aggr looking at tcga and local data (except the modelling data)
sign_node_d_all, sign_node_d_all_data = significant_loci_all_data(node_d_aggr_tcga, node_d_aggr, node_d_nonaggr_tcga,
                                                                  node_d_nonaggr, features_names_node_d)
sign_btw_c_all, sign_btw_c_all_data = significant_loci_all_data(btw_c_aggr_tcga, btw_c_aggr, btw_c_nonaggr_tcga,
                                                                btw_c_nonaggr, features_names_btw)
sign_dg_c_all, sign_dg_c_all_data = significant_loci_all_data(dg_c_aggr_tcga, dg_c_aggr, dg_c_nonaggr_tcga,
                                                              dg_c_nonaggr, features_names_dg)

node_d_all_melt = pd.concat([node_d_melt, node_d_melt_tcga])
sign_node_d_all_melt = node_d_all_melt[node_d_all_melt["variable"].isin(sign_node_d_all)]

btw_c_all_melt = pd.concat([btw_c_melt, btw_c_melt_tcga])
sign_btw_c_all_melt = btw_c_all_melt[btw_c_all_melt["variable"].isin(sign_btw_c_all)]

dg_c_all_melt = pd.concat([dg_c_melt, dg_c_melt_tcga])
sign_dg_c_all_melt = dg_c_all_melt[dg_c_all_melt["variable"].isin(sign_dg_c_all)]

sns.set(font_scale=1.5)

plt.figure(figsize=(15, 20))
sns.boxplot(data=sign_node_d_all_melt, x="value", y="variable", hue="aggressive", ).set_title(
    "node degree\nall samples (AE and TCGA)\nloci with significant difference")
# plt.savefig("results/pnet_plots/" + base + "/boxplot_node_d_all_data_sign" + "_" + str(n) + "_" + str(t) + "_" + str(k) + ".png")
# plt.savefig("results/pnet_plots/" + base + "/boxplot_node_d_all_data_sign" + "_" + str(n) + "_" + str(t) + "_" + str(k) + ".pdf")
plt.show(block=True)
# plt.close()

plt.figure(figsize=(15, 20))
sns.boxplot(data=sign_btw_c_all_melt, x="value", y="variable", hue="aggressive").set_title(
    "betweenness centrality\nall samples (AE and TCGA)\nloci with significant difference")
# plt.savefig("results/pnet_plots/" + base + "/boxplot_btw_c_all_data_sign" + "_" + str(n) + "_" + str(t) + "_" + str(k) + ".png")
# plt.savefig("results/pnet_plots/" + base + "/boxplot_btw_c_all_data_sign" + "_" + str(n) + "_" + str(t) + "_" + str(k) + ".pdf")
plt.show(block=True)
# plt.close()

plt.figure(figsize=(15, 20))
sns.boxplot(data=sign_dg_c_all_melt, x="value", y="variable", hue="aggressive").set_title(
    "degree centrality\nall samples (AE and TCGA)\nloci with significant difference")
# plt.savefig("results/pnet_plots/" + base + "/boxplot_dg_c_all_data_sign" + "_" + str(n) + "_" + str(t) + "_" + str(k) + ".png")
# plt.savefig("results/pnet_plots/" + base + "/boxplot_dg_c_all_data_sign" + "_" + str(n) + "_" + str(t) + "_" + str(k) + ".pdf")
plt.show(block=True)
# plt.close()
