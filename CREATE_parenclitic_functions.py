"""
This script includes the functions necessary to create the parenclitic networks.
These functions are used in parenclitic.py
"""

import numpy as np
import pandas as pd
import networkx as nx
from sklearn.mixture import BayesianGaussianMixture
from scipy.spatial.distance import mahalanobis

def coord_matrix(beta):
    # transform the beta data in a coordinates array form to be used for linear regression and distance calculation
    # 4 dimensional locus i, locus j, sample s, coordinates(x,y)
    # the diagonal (i,i) and below (j,i) with j>i will not be assigned a value (remain 0) due to symmetry of the array
    # input: df lxk with beta values
    # output: coordinates array for each pair of loci for each sample, naming of the loci and samples in the array

    l, k = beta.shape
    
    # keep record of the locus name per index - name should be the same for any beta input (aggr and non-aggr)
    locus_naming = [beta.index[i] for i in range(l)]
    sample_naming = [beta.columns[i] for i in range(k)]

    # initialise coordinates array
    b_coord = np.zeros((l, l, k, 2))

    for i in range(l):
      for j in range(i + 1, l):
          for s in range(k):
              b_coord[i, j, s, 0] = beta[sample_naming[s]][locus_naming[i]]
              b_coord[i, j, s, 1] = beta[sample_naming[s]][locus_naming[j]]

    return b_coord, locus_naming, sample_naming


def calculate_distance(b_coord_ref, b_coord, l, k, locus_naming):
    # this section calculates the edge weight from the reference for each subject separately
    # i.e. distance of an observation from the reference model
    # initialise weight array w; lxlxk -> the diagonal and lower half are 0 (symmetry - only upper half calculated)
    # y_pred[i,j,s] will store the predicted value for locus j with from the input locus i for the sample s
    # input: beta coordinates for the reference data and the class to be compared to the ref model
    # l: the number of loci, k: the number of samples
    # output: array w lxlxk, w[i,j,k]: the weight between loci i and j of sample s

    # initialise weight and predicted values arrays
    w = np.zeros((l, l, k))
    y_pred = np.zeros((l, l, k))
    for i in range(l):
        for j in range(i + 1, l):
            model = BayesianGaussianMixture(n_components=4, random_state=10).fit(b_coord_ref[i, j, :, :])

            for s in range(k):
                # predict the beta value at locus j using the reference model with input the beta value at locus i
                # X = b_coord[i, j, s, 0]
                # y = b_coord[i, j, s, 1]
                y_pred[i, j, s] = model.predict(b_coord[i, j, s, :].reshape(1, -1))

                # predict to which gaussian the sample belongs
                c = y_pred[i,j,s].astype("int")
                # compute the inverse covariant matrix
                icv = np.linalg.inv(model.covariances_[c])
                # compute the mahalanobis distance
                w[i, j, s] = mahalanobis(b_coord[i,j,s,:], model.means_[c], icv)

    return w


def weights_list(w, threshold, locus_naming_target, locus_naming_ref):
    # creates a list of tuples with loci pairs and corresponding weight
    # this list can be used to create the graph
    # input: weights array, threshold for edge weight, locus naming of both classes for sanity check (data correctness)
    # ref corresponds to the data used for the reference model
    # target corresponds to the class for which we want to calculate the parenclisis from the reference model
    # output: list of tuples (locus name i, locus name j, weight of edge connecting loci i and j)

    # initialise list
    g_summary = []
    l, l, k = w.shape

    # check locus_naming_target == locus_naming_ref first
    # the names must coincide
    # w doesn't include the naming only the indexing
    if locus_naming_target == locus_naming_ref:
        print("locus namings coincide: the list can be created")

        # create list of weights
        # the graph will be undirected hence direction (i->j or j->i) doesn't play a role
        for s in range(k):
            weights_list = []
            for i in range(l):
                for j in range(i + 1, l):
                    # only include non 0 weights above the diagonal
                    if w[i, j, s] > threshold:
                        weights_list = weights_list + [(locus_naming_target[i], locus_naming_target[j], w[i, j, s])]
            g_summary.append(weights_list)
    else:
        print(" locus names don't coincide: check the naming and indexing before proceeding")

    return g_summary


def create_graph(k, locus_naming, g_summary, base, title, name, l, t):
    # creates and saves graph, and prints summary statistics
    # k = k_aggr or k_nonaggr, locus_naming = locus_naming_aggr or locus_naming_non_aggr
    # g_summary = g_aggr_summary or g_nonaggr_summary, base = "based_on_tcga", title = "aggressive" or "non aggressive"
    # name = aggr or non_aggr # sample name, l number of loci, t is the threshold
    # saves topological features

    # initialise dataframes for topo measures
    node_d = pd.DataFrame(data=None, index=name, columns=locus_naming, dtype=None, copy=False)
    d_c = pd.DataFrame(data=None, index=name, columns=locus_naming, dtype=None, copy=False)
    btw_c = pd.DataFrame(data=None, index=name, columns=locus_naming, dtype=None, copy=False)
    
    # create graphs
    # aggressive
    # individual plots
    for s in range(k):
        # create graph
        G = nx.Graph()
        G.add_nodes_from(locus_naming)
        G.add_weighted_edges_from(g_summary[s])
        
        # save graph
        # nx.write_gml(G, "~/Results/pnets/" + base + "/graphs/pnet_" + str(l) + "_" + str(t) + "_"
        #              + name[s] + "_" + title + ".gml")
        nx.write_graphml_lxml(G, "~/Results/pnets/" + base + "/graphs/pnet_" + str(l) + "_" + str(t) +
                              "_" + name[s] + "_" + title + ".graphml")

        # summary statistics
        print("--------------")
        print("Stats for sample " + name[s])
        print("Degree")
        stat = pd.DataFrame(G.degree)
        print(stat.describe())
        print("number of nodes with 0 degree: " + str(sum(stat[1] == 0)))
        print("Edges")
        print(pd.DataFrame(G.edges).describe())
        print("--------------")

        # computing topological features
        node_degrees = dict(G.degree())
        node_d.loc[name[s]][node_degrees.keys()] = list(node_degrees.values())
        d_centrality = dict(nx.degree_centrality(G))
        d_c.loc[name[s]][d_centrality.keys()] = list(d_centrality.values())
        btw_centr = nx.betweenness_centrality(G, normalized=True)
        btw_c.loc[name[s]][btw_centr.keys()] = list(btw_centr.values())

    node_d.to_csv("~/Data/" + base + "/node_degrees_" + title + "_" + str(l) + "_" + str(t) + ".csv")
    d_c.to_csv("~/Data/" + base + "/degree_centrality_" + title + "_" + str(l) + "_" + str(t) + ".csv")
    btw_c.to_csv("~/Data/" + base + "/betweenness_centrality_" + title + "_" + str(l) + "_" + str(t) + ".csv")




