"""
This script creates the baseline model for non-aggressive samples and
then for each sample not used in the modelling creates the parenclitic network.
"""

import os
import pandas as pd

from parenclitic_functions import coord_matrix, calculate_distance, weights_list, create_graph


# load beta values and sample information
print("loading data")
# point to the right directory
os.chdir("~/Data")
# my data
beta = pd.read_csv("beta_bmiq_experimental.csv", index_col=0)
targets = pd.read_csv("Methylation/all_idat/sample_sheet.csv", skiprows=7, index_col="barcode")
#TCGA
beta_tcga = pd.read_csv("Methylation/TCGA/TCGA_beta_bmiq_experimental_renamed.csv", index_col=0)
targets_tcga = pd.read_excel("mmc3.xls", skiprows=2, index_col="Sample ID")
targets_tcga.rename(columns={"Clinically Aggressive (extensive invasion, "
                             "local regional disease, local recurrence)" : "Clinically_Aggressive"}, inplace=True)


# 850k and 450k overlap
overlap = [gene for gene in beta_tcga.index if gene in beta.index]
print("overlap of topvar genes in my data and TCGA: " + str(len(overlap)))
beta = beta.loc[overlap,]
beta_tcga = beta_tcga.loc[overlap,]

# # compute variance per site and sort for the TCGA samples
# # choose top l most variant genes
print("computing variance per site")
l = 1000
var = beta_tcga.var(axis=1)
var = var.sort_values(ascending=False)
topvar = var[0:l]
beta_topvar_tcga = beta_tcga.loc[topvar.index,:]
beta_tcga = beta_topvar_tcga
beta = beta.loc[topvar.index,:]

# create separate dataframes for aggressive and non-aggressive
# and choose the topvar sub-dataset
print(" creating aggr & non_aggr dataframes")
non_aggr = targets.index[targets.Clinically_Aggressive == "no"]
aggr = targets.index[targets.Clinically_Aggressive == "yes"]
non_aggr_tcga = targets_tcga.index[targets_tcga.Clinically_Aggressive == "no"]
aggr_tcga = targets_tcga.index[targets_tcga.Clinically_Aggressive == "yes"]

# keep 16 random non_aggr tcga samples for training and use the rest for the model
random_samples_index = [10, 85, 86, 93, 23, 154, 52, 44, 136, 138, 112, 101, 140, 151, 33, 35]
model_samples_index = [i for i in range(len(non_aggr_tcga)) if i not in random_samples_index]

b_nonaggr = beta[non_aggr.astype("string")]
b_aggr = beta[aggr.astype("string")]
b_nonaggr_tcga_model = beta_tcga[non_aggr_tcga.astype("string")[model_samples_index]]
b_nonaggr_tcga_train = beta_tcga[non_aggr_tcga.astype("string")[random_samples_index]]
b_aggr_tcga = beta_tcga[aggr_tcga.astype("string")]
print(b_aggr.shape)

# transform the data into coordinates array to be used for linear regression and distance calculation
# 4 dimensional locus i, locus j, sample s, coordinates(x,y)
# the diagonal (i,i) and below (j,i) with j>i will not be assigned a value (remain 0) due to the symmetry of the array
print("calling coordinates function")

# Array Express
# non-aggressive
b_nonaggr_coord, gene_naming_non_aggr, sample_naming_non_aggr = coord_matrix(b_nonaggr)
# aggressive
b_aggr_coord, gene_naming_aggr, sample_naming_aggr = coord_matrix(b_aggr)
print("array express aggr and non aggr coordinates created")

# TCGA
# non-aggressive model
b_nonaggr_coord_tcga_model, gene_naming_non_aggr_tcga_model, sample_naming_non_aggr_tcga_model = coord_matrix(
    b_nonaggr_tcga_model)
# non-aggressive train
b_nonaggr_coord_tcga_train, gene_naming_non_aggr_tcga_train, sample_naming_non_aggr_tcga_train = coord_matrix(
    b_nonaggr_tcga_train)
# aggressive
b_aggr_coord_tcga, gene_naming_aggr_tcga, sample_naming_aggr_tcga = coord_matrix(b_aggr_tcga)
print("tcga aggr and non aggr coordinates created")

k_aggr = b_aggr.shape[1]
k_nonaggr = b_nonaggr.shape[1]
k_aggr_tcga = b_aggr_tcga.shape[1]
k_nonaggr_tcga_model = b_nonaggr_tcga_model.shape[1]
k_nonaggr_tcga_train = b_nonaggr_tcga_train.shape[1]

# calculate the edge weight, i.e. distance of a sample from the expected non-aggressive model
print("calculating weights for aggr and non_aggr")

base = "based_on_tcga_data"
w_aggr = calculate_distance(b_nonaggr_coord_tcga_model, b_aggr_coord, l, k_aggr, gene_naming_aggr)
w_nonaggr = calculate_distance(b_nonaggr_coord_tcga_model, b_nonaggr_coord, l, k_nonaggr, gene_naming_non_aggr)
w_aggr_tcga = calculate_distance(b_nonaggr_coord_tcga_model, b_aggr_coord_tcga, l, k_aggr_tcga, gene_naming_aggr_tcga)
w_nonaggr_tcga = calculate_distance(b_nonaggr_coord_tcga_model, b_nonaggr_coord_tcga_train, l, k_nonaggr_tcga_train,
                                    gene_naming_non_aggr_tcga_train)
print("Distance calculated ")

# create a list of tuple with loci pairs and corresponding weight
thresh = 2
g_nonaggr_summary = weights_list(w_nonaggr, thresh, gene_naming_non_aggr, gene_naming_non_aggr_tcga_model)
g_aggr_summary = weights_list(w_aggr, thresh, gene_naming_aggr, gene_naming_non_aggr_tcga_model)
g_nonaggr_summary_tcga = weights_list(w_nonaggr_tcga, thresh, gene_naming_non_aggr_tcga_train,
                                      gene_naming_non_aggr_tcga_model)
g_aggr_summary_tcga = weights_list(w_aggr_tcga, thresh, gene_naming_aggr_tcga, gene_naming_non_aggr_tcga_model)

print("weights in correct format - creating graphs")

# create graphs
# Array Express
create_graph(k_aggr, gene_naming_aggr, g_aggr_summary, base, "aggressive", aggr, l, t)
create_graph(k_nonaggr, gene_naming_non_aggr, g_nonaggr_summary, base, "non_aggressive", non_aggr, l, t)
# TCGA
create_graph(k_aggr_tcga, gene_naming_aggr_tcga, g_aggr_summary_tcga, base, "aggressive_tcga", aggr_tcga, l, t)
create_graph(k_nonaggr_tcga_train, gene_naming_non_aggr_tcga_train, g_nonaggr_summary_tcga, base, "non_aggressive_tcga",
             non_aggr_tcga[random_samples_index], l, t)





