"""
This script reads the graph files (.graphml) of all samples and calculates their topological features.
The topological features are saved in csv files.
"""

import os
import networkx as nx
import pandas as pd
import time

n = 1000
k = 2
base = "based_on_tcga_data"

graph_files = os.listdir("~/Results/pnets/" + base + "/graphs/")
graph_files = [g for g in graph_files if str(n) + "_" + str(k) in g]

print(graph_files)

# initialise graph and features dictionary
influence = {}  # influential nodes
node_d = {}  # node degree
eig_c = {}  # eigenvector centrality
btw_c = {}  # betweenness centrality
dg_c = {}   # degree centrality
sec_ord_c = {}  # second order centrality
ecc = {}    # eccentricity

number_of_features = 7
number_of_samples = len(graph_files)

for f in graph_files:
    sample_name = sample_name = f[11:-8]
    print(sample_name)

    start_time = time.perf_counter()
    G = nx.read_graphml("~/Results/pnets/" + base + "/graphs/" + f)

    influence[sample_name] = nx.voterank(G, number_of_nodes=200)
    node_d[sample_name] = pd.DataFrame.from_dict(G.degree()).set_index(0)
    eig_c[sample_name] = nx.eigenvector_centrality_numpy(G, max_iter=20)
    btw_c[sample_name] = nx.betweenness_centrality(G, normalized=True)
    dg_c[sample_name] = nx.degree_centrality(G)
    
    # topo features for largest subgraph
    G_sub = sorted(nx.connected_components(G), key=len, reverse=True)
    G_sub0 = G.subgraph(G_sub[0])
    sec_ord_c[sample_name] = nx.second_order_centrality(G_sub0)
    ecc[sample_name] = nx.eccentricity(G_sub0)

    finish_time = time.perf_counter()
    print("sample " + sample_name + " and features done in " + str(finish_time-start_time))
    

# dict to df
print("converting to dataframes")
node_d_df = pd.concat(node_d, axis=1)
node_d_df.columns = node_d_df.columns._get_level_values(0)
eig_c_df = pd.DataFrame.from_dict(eig_c)
btw_c_df = pd.DataFrame.from_dict(btw_c)
dg_c_df = pd.DataFrame.from_dict(dg_c)
sec_ord_c_df = pd.DataFrame.from_dict(sec_ord_c)
ecc_df = pd.DataFrame.from_dict(ecc)

# save topo features (tf)
print("saving to csv files")
node_d_df.to_csv("~/Data/" + base + "/tf_node_degrees_all_" + str(n) + "_" + str(k) + ".csv")
eig_c_df.to_csv("~/Data/" + base + "/tf_eigenvector_centrality_all_" + str(n) + "_" + str(k) + ".csv")
btw_c_df.to_csv("~/Data/" + base + "/tf_betweenness_centrality_all_" + str(n) + "_" + str(k) + ".csv")
dg_c_df.to_csv("~/Data/" + base + "/tf_degree_centrality_all_" + str(n) + "_" + str(k) + ".csv")
sec_ord_c_df.to_csv("~/Data/" + base + "/tf_second_order_centrality_all_" + str(n) + "_" + str(k) + ".csv")
ecc_df.to_csv("~/Data/" + base + "/tf_eccentricity_all_" + str(n) + "_" + str(k) + ".csv")






