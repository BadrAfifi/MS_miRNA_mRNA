BiocManager::install("apeglm")
install.packages("devtools")
library(devtools)
install_github("stephens999/ashr")
devtools::install_github('colearendt/xlsx')

BiocManager::install("org.Hs.eg.db")

library(xlsx)
library(org.Hs.eg.db)
library(biomaRt)
library(ashr)
library(apeglm)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(ashr)

bico
gen_counts <- read.table("GSE172009_Raw_gene_counts_matrix.txt", header = T, row.names =  1)
str(gen_counts)

# Create condition vector
Condition <- c(rep("HC", 4),rep("RRMS", 4))

# Create data frame
metadata <- data.frame(Condition )

# Assign the row names of the data frame
rownames(metadata) <- colnames(gen_counts)
all(rownames(metadata) == colnames(gen_counts))

# If "FULSE":
idx <- match(colnames(gen_counts), rownames(metadata))
ordered_metaD <- metadata[idx, ]

all(rownames(ordered_metaD) == colnames(gen_counts))

# then 
#           --- DESeq2 ----

#Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData= round(gen_counts), 
                              colData= metadata, 
                              design= ~ Condition) # tidy = F


#Variations to the standard workflow
#Normalization Counts : calculation
dds <- estimateSizeFactors(dds)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

sizeFactors(dds)

normalized_counts <- counts(dds, normalized = T)

#Unsupervised clustering analyses: log transformation
vsd_dds <- vst(dds, blind=TRUE)

#  1 ----- Hierarchical clustering with correlation heatmaps

# Extract the vst matrix from the object
vsd_mat <- assay(vsd_dds)

# Compute pairwise correlation values
vsd_cor <- cor(vsd_mat)

# Plot heatmap
pheatmap(vsd_cor, annotation = select(metadata, Condition))

#  2 ----- Principal component analysis (PCA):


# Plot PCA
plotPCA(vsd_dds, intgroup="Condition")






# ----- DESeq2workflow:Running;

dds2 <- DESeq(dds)


#DESeq2model-mean-variancerelationship:

#Calculating mean for each gene (each row)
mean_counts <- apply(gen_counts[,1:3],1, mean)

# Calculating variance for each gene (each row)
variance_counts <- apply(gen_counts[,1:3],1, var)

# Creating data frame with mean and variance for every gene
df <- data.frame(mean_counts, variance_counts)

ggplot(df)+geom_point(aes(x=mean_counts, y=variance_counts))+
  scale_y_log10()+scale_x_log10()+
  xlab("Mean counts per gene")+ylab("Variance per gene")

# ---- Plot dispersion estimates
plotDispEsts(dds2)


# --- DESeq2 contrasts
resultsNames(dds2)
res <- results(dds2,
               contrast = c("Condition", "RRMS","HC"),
               alpha = 0.05,lfcThreshold = 2,
               cooksCutoff= TRUE) 
summary(res)
mcols(res)
#contrast = c("condition_factor","level_to_compare","base_level")
#cooksCutoff; T is deafult.
plotMA(res, ylim=c(-8,8))


# ---- DESeq2 LFC shrinkage estimators

# --- 1- ashr : using BH  produces p adj 
resAshT <- lfcShrink(dds2 ,
                 contrast = c("Condition", "RRMS","HC"),
                 res= res,
                 coef=2,
                 type="ashr")
mcols(resAshT)
head(resAshT, n=10)
summary(resAshT)

plotMA(resAshT, ylim=c(-3,3), cex=.8)
abline(h=c(-1,1), col="dodgerblue", lwd=2)


# OR
# --- 2- apeglm : produces s-value
resApeT <- lfcShrink(dds2 ,
                     lfcThreshold = 1,
                     res= res,
                     coef=2,
                     type="apeglm")

plotMA(resApeT, ylim=c(-3,3), cex=.8)
abline(h=c(-1,1), col="green", lwd=2)

mcols(resApeT)
head(resApeT, n=10)
summary(resApeT)


# Finally :

sig_rse <- subset(resAshT, padj < 0.05 & abs(log2FoldChange) > 1)
sig_rse <- data.frame(sig_rse[order(sig_rse$padj),]) # ordering

# Getting gene symbols for Transcripts IDs
gene_Refseq_IDs <- gsub("\\..*", "", rownames(sig_rse))
keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)
cols <- c("SYMBOL", "GENENAME")
ensids <- gene_Refseq_IDs
gene_info <- select(org.Hs.eg.db, keys=ensids, columns=cols, keytype="REFSEQ")


sig_rse <- bind_cols(rownames(sig_rse),sig_rse,gene_info)
names(sig_rse)[1] = "Transcript ID"
sig_rse_up <- subset(sig_rse,log2FoldChange > 1)
sig_rse_down <- subset(sig_rse,log2FoldChange < 1)


write.table(bind_cols(sig_rse$`Transcript ID`,sig_rse$SYMBOL,sig_rse$GENENAME),
            file = "DEGs.txt",row.names = F,col.names = F,quote = F) 
write.table(bind_cols(sig_rse_up$SYMBOL),
            file = "DEGs_Up.txt",row.names = F,col.names = F,quote = F) 
write.table(bind_cols(sig_rse_down$SYMBOL),
            file = "DEGs_down.txt",row.names = F,col.names = F,quote = F)

# collecting NAN symbols
sig_rse <- sig_rse[order(sig_rse$SYMBOL),]
NAN_genes <- gsub("\\..*", "", rownames(sig_rse))
write.table(bind_cols(sig_rse$`Transcript ID`,sig_rse$SYMBOL),
            file = "NAN_DEGs.txt",row.names = F,col.names = F,quote = F) 

