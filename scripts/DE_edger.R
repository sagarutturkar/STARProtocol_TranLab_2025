library(DESeq2)
library(ggplot2)
library(gplots)
library(tidyverse)
library(RColorBrewer)
library(edgeR)
library(ggrepel)
library(ComplexHeatmap)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

countfile     = args[1] #path for the counts file
control       = args[2] #name for control (as in counts file)
treatment     = args[3] #name for treatment (as in counts file)
control_rep   = args[4] #number of replicates for control
treatment_rep = args[5] #number of replicates for treatment
path          = args[6] #path for annotation file

outprefix <- paste(treatment, control, sep="_vs_")
dir.create(outprefix)
setwd(outprefix)

header = paste(rep("#", 50), collapse = "")

sink(file = paste0(outprefix,".sessioninfo.txt"))

cat(paste(header, "#Version Information", header, sep = '\n'))
cat('\n')
version
cat('\n')

cat(paste(header, "#Session Information", header, sep = '\n'))
cat('\n')
sessionInfo()
sink()

#Save log to file

sink(file = paste0(outprefix,".log.txt"))

Anno <- read.table(path, sep = "\t", stringsAsFactors = F, header = TRUE, quote = "")

raw.data = read.table(countfile, row.names="Gene_ID", header = T, 
                      sep = "\t", as.is = T, check.names = F, quote = "")

control
treatment

group <- factor(rep(c(control, treatment),
                    times=c(control_rep, treatment_rep)),
                levels = c(control, treatment))

group
control_mat<- select(raw.data,starts_with(control))
treatment_mat<- select(raw.data,starts_with(treatment))
raw.counts<- merge(control_mat,treatment_mat,by="row.names",all.x=TRUE)

raw.counts<- column_to_rownames(raw.counts, var="Row.names")


###################################################################################
#                          edgeR Analysis
###################################################################################

y <- DGEList(counts=as.matrix(raw.counts), group=group)
x <- calcNormFactors(y)


keep <- filterByExpr(y)
filtered.data <- y[keep,keep.lib.sizes=FALSE]
y <- calcNormFactors(filtered.data)


design <- model.matrix(~0+group)


colnames(design)<- levels(y$samples$group)
y <- estimateDisp(filtered.data,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,contrast=c(-1,1))

edgeR_DE = as.data.frame(topTags(qlf, sort.by = "PValue", n = Inf))
edgeR_DE = rownames_to_column(edgeR_DE, "Gene_ID")
edgeR_DE = left_join(edgeR_DE, Anno, by="Gene_ID")

write.table(edgeR_DE,        file = paste0(outprefix,".DE_edgeR_All.txt"), sep = "\t", quote = F, row.names = F)

saveRDS(edgeR_DE, file = paste0(outprefix,"edgeR_DE.rds"))


###################################################################################
#                          Filtering
###################################################################################

# EdgeR - FDR 5% filtered
edgeR_DE_FDR5P = dplyr::filter(edgeR_DE, FDR < 0.05)
dim(edgeR_DE_FDR5P)

write.table(edgeR_DE_FDR5P,  file = paste0(outprefix,".DE_edgeR_5PCT.txt"), sep = "\t", quote = F, row.names = F)

sink()
