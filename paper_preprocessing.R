# This script includes the preprocessing and QC of the idat files
# and results in the final beta values used for further analysis.

library(minfi)
library(wateRmelon)

# point to the right directory
setwd('~/Data/Methylation/')

# point targets to the folder with the idat files and the samplesheet
targets <- read.metharray.sheet( file.path( './all_idat'), verbose= TRUE)
RGset <- read.metharray.exp(targets = targets, verbose=TRUE, force=TRUE)
anno = getAnnotation(RGset)


# samples with median methylated and unmethylated intensities <9.5 to be removed
col_median_green = colMedians(as.matrix(RGset@assays@data@listData$Green))
col_median_red = colMedians(as.matrix(RGset@assays@data@listData$Red))
sum(col_median_red <9.5)    # 0 -> no samples to be removed
sum(col_median_green <9.5)  # 0 -> no samples to be removed

remove(col_median_green, col_median_red)

# probes with detection p-value>0.05 regarded as failed
detp = detectionP(RGset)
table(detp>0.05)

# samples with >10% failed probes to be removed
# non such samples -> no samples to be removed
failed_samples = c()
for (i in 1:ncol(RGset@assays@data@listData$Green)) {
  n_failed = sum(RGset@assays@data@listData$Green[,i]<9.5)
  n = length(RGset@assays@data@listData$Green[,i])
  ratio_failed = n_failed/n
  # print(n_failed)
  if (ratio_failed>0.1) {
    print(c(i," failed"))
    failed_samples =c(failed_samples, RGset@assays@data@listData$Green[,i])
    
  }
} 

# failure considered based on pval>0.05
# failed_probes list of names of probes to be set to NA 
failed_probes = c()
for (i in 1:nrow(detp)) {
  # print(i)
  n_failed = sum(detp[i,]>0.05)
  n = length(detp[i,])
  ratio_failed = n_failed/n
  if (ratio_failed>0.1) {
    failed_probes =c(failed_probes,rownames(detp)[i])
  }
} 

remove(detp, i, n, n_failed, ratio_failed)

# transform to Mset
# needs to be converted to Mset to use the same row names as detp

# background intensity correction & dye bias correction
Mset = preprocessNoob(RGset, offset = 15, dyeCorr = TRUE, verbose = FALSE,
                      dyeMethod="single")
# beta values after probe bias correction using BMIQ
b_bmiq= BMIQ(Mset)

# remove failed probes
Mset = Mset[!(rownames(Mset) %in% failed_probes)]

remove(RGset, failed_samples, failed_probes)
gc()

# drop SNP probes (probes that measure SNPs, not probes containing SNPs) 
# and CH probes (non-CpG methylation)
Mset = dropMethylationLoci(Mset, dropRS = TRUE, dropCH = TRUE)

# remove SNP-related probes
# it needs to be converted to GMset first to use dropLociWithSnps
GMset = mapToGenome(Mset)
GMset = dropLociWithSnps(GMset)

remove(Mset)
# remove sex probes
non_sex_probes =  rownames(anno)[!(anno@listData$chr %in% c("chrX", "chrY"))]
GMset = GMset[rownames(GMset) %in% non_sex_probes,]

remove(non_sex_probes)

# remove the SNP, non-CpG and sex probes from the BMIQ normalised beta values
b_bmiq = b_bmiq[rownames(b_bmiq) %in% rownames(GMset),]
# write.csv(b_bmiq,file="beta_bmiq_experimental.csv")


