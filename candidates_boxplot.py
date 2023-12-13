"""
This script generates the boxplots comparing the beta values between aggressive and non-aggressive tumours
at the candidate loci.
"""

import os
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
from statannot import add_stat_annotation

# point to the right directory
os.chdir('~/Data/Methylation/')

b = pd.read_csv("./beta_bmiq_experimental.csv", index_col=0)

samples = pd.read_csv("./all_idat/sample_sheet.csv", index_col="barcode", header=7)
ClA = samples.Clinically_Aggressive.dropna()

candidates = ["cg10928544", "cg22933800", "cg24435747", "cg05146756", "cg03758477",
              "cg13563298", "cg22891070", "cg22332722", "cg26950867", "cg27576485"]

b_candidates = b.loc[candidates, ClA.index]
samples = samples.loc[ClA.index]
b_candidates_and_phenotype = pd.concat([b_candidates.transpose(), ClA], axis=1, join="inner")

# define the axes and the layout of the figure
a = [(0,0), (0,1), (1,0), (1,1), (2,0), (2,1), (3,0), (3,1), (4,0), (4,1),]
i = 0
f, axes = plt.subplots(5, 2, layout="constrained")
f.set_figheight(15)
f.set_figwidth(8)
f.suptitle("Beta values comparison in candidate CG loci (Mann-Whitney test)")

for cpg in candidates:
    ax = sns.boxplot(data=b_candidates_and_phenotype, y=cpg, x="Clinically_Aggressive", ax = axes[a[i]])
    test_results = add_stat_annotation(ax, plot="boxplot", data=b_candidates_and_phenotype, y=cpg, x="Clinically_Aggressive",
                                       box_pairs=[("yes", "no")],
                                       test="Mann-Whitney",
                                       text_format='star',
                                       loc='inside', verbose=2)
    i += 1

plt.savefig("./results/pnet_plots/candidates_boxplots.png")
plt.savefig("./results/pnet_plots/candidates_boxplots.pdf")