# This script performs differential methylation analysis 
# using DMPfinder and bumphunter 
# both between aggressive/non-aggressive and metastatic/non-metastatic
# and compares the results in a venn diagramm.

library(minfi)
library(ggvenn)
library(bumphunter)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# point to the right directory
setwd('~/Data/Methylation/')

# read data
b = read.csv("~/Data/Methylation/beta_bmiq_experimental.csv",
             check.names = FALSE)
rownames(b) = b[,1]
b = b[,2:54]
targets <- read.metharray.sheet( file.path( './all_idat'), verbose= TRUE)
anno = read.csv("~/Data/Methylation/anno_minfi.csv")
rownames(anno) = anno$X
anno= anno[,2:length(colnames(anno))]

# Data frame with column annotations.
mat = data.frame(group = targets)

# define attribute of interest for the analysis
Group = mat$group.Sample_Group
Metastasis = mat$group.Metastatic
ClAggr = mat$group.Clinically_Aggressive


#mat_col = data.frame(interest)
mat_col = data.frame(Group = Group, Met = Metastasis, ClA = ClAggr,
                     cluster=c(1:length(colnames(b)))*0)
rownames(mat_col) <- mat$group.barcode

# add the cluster to the dataframe according to the group 
for (i in c(1:53)){
  if (mat_col$Group[i] %in% c("HRAS", "RET", "NF1")) { mat_col$cluster[i] = "C2"}
  else if (mat_col$Group[i] %in% c("SDHA", "SDHB", "SDHD", "IDH3B")) {
    mat_col$cluster[i] = "C1A"}
  else if (mat_col$Group[i] %in% c("VHL")) { mat_col$cluster[i] = "C1B"}
  else if (mat_col$Group[i] %in% c("normal")) { mat_col$cluster[i] = "normal"}
  else if (mat_col$Group[i] %in% c("no_mut")) { mat_col$cluster[i] = "sporadic"}
}

###############################################################################

# phenotypes of interest
met = mat_col$Met[mat_col$Met %in% c("yes", "no")]
aggr = mat_col$ClA[mat_col$ClA %in% c("yes", "no")]

# SDHB.met = (mat_col$Group == "SDHB") & (mat_col$Met %in% c("yes", "no"))
# SDHB.cla = (mat_col$Group == "SDHB") & (mat_col$ClA %in% c("yes", "no"))

b = as.matrix(b)

# differentially methylated positions
dmp.met = dmpFinder(b[,mat_col$Met %in% c("yes", "no")], pheno= met,
                    type = 'categorical')
dmp.aggr = dmpFinder(b[,mat_col$ClA %in% c("yes", "no")], pheno= aggr,
                     type = 'categorical')


significant_dmp = function(dmp, q, p){
  dmp.lowq = dmp[dmp$qval<q, ]
  dmp.lowp = dmp.lowq[dmp.lowq$pval<p,]
  
  return(dmp.lowp)
  
}

par(mfrow = c(5,5))
dmsites.met = significant_dmp(dmp = dmp.met, q = 0.05, p = 0.005)
dmsites.cla = significant_dmp(dmp = dmp.aggr, q = 0.05, p = 0.005)

# get the gene names for the cg loci
dmp.genes.met = anno$UCSC_RefGene_Name[rownames(anno) %in% rownames(dmsites.met)]
dmp.genes.met = gsub("(;).*","",dmp.genes.met)

dmp.genes.cla = anno$UCSC_RefGene_Name[rownames(anno) %in% rownames(dmsites.cla)]
dmp.genes.cla = gsub("(;).*","",dmp.genes.cla)

b.met = b[,mat_col$Met %in% c("yes", "no")]
b.clA = b[,mat_col$ClA %in% c("yes", "no")]

attr.met = mat_col$Met[mat_col$Met %in% c("yes", "no")]
attr.clA = mat_col$Met[mat_col$ClA %in% c("yes", "no")]

###############################################################################
# differentially methylated regions
met.mat = model.matrix(~ met)
aggr.mat = model.matrix(~ aggr)

dmrs.met = bumphunter(as.matrix(b.met),
                      design = met.mat, chr = anno$chr, pos = anno$pos,
                      pickCutoff = TRUE, # cutoff = 0.3,
                      B=100, type = 'Beta',
                      nullMethod = 'bootstrap', verbose = TRUE)
dmrs.met$table = dmrs.met$table[dmrs.met$table$p.value<0.05,]

dmrs.aggr = bumphunter(as.matrix(b.clA),
                       design = aggr.mat, chr = anno$chr, pos = anno$pos,
                       pickCutoff = TRUE, # cutoff = 0.3,
                       B=100, type = 'Beta',
                       nullMethod = 'bootstrap', verbose = TRUE)
dmrs.aggr$table = dmrs.aggr$table[dmrs.aggr$table$p.value<0.05,]
gc()

# identify genes detected with bumphunter
transcript.anno = annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                      annotation="org.Hs.eg.db")
bump.genes.met = matchGenes(dmrs.met$table, transcript.anno)
bump.genes.aggr = matchGenes(dmrs.aggr$table, transcript.anno)

# identify dm sites detected with bumphunter
cpg.met = rownames(anno)[anno$pos %in% dmrs.met$table$start]
cpg.met = cpg.met[!is.na(cpg.met)]
# write.csv(cpg.met, file="./results/DM/DMsites_bumphunter_metastatic.csv")

cpg.aggr = rownames(anno)[anno$pos %in% dmrs.aggr$table$start]
cpg.aggr = cpg.aggr[!is.na(cpg.aggr)]
# write.csv(cpg.aggr, file="./results/DM/DMsites_bumphunter_aggressive.csv")

b.bump.met = b[rownames(b) %in% cpg.met,]
b.bump.aggr = b[rownames(b) %in% cpg.aggr,]

par(mfrow = c(1,1))
# png(filename = "./results/DM/DMP_all_venn.png",
#     width = 900, height = 600, units = "px", pointsize = 12)
ggvenn(list(metastatic.dmpFinder = unique(dmp.genes.met),
                       aggressive.dmpFinder = unique(dmp.genes.cla),
                       metastatic.bumphunter = unique(bump.genes.met$name),
                       aggressive.bumphunter = unique(bump.genes.aggr$name)),
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
       stroke_size = 0.5,
       set_name_size = 6,
       text_size =5,
       show_percentage = FALSE) + scale_x_continuous(expand = expansion(mult = .2))
# dev.off()

