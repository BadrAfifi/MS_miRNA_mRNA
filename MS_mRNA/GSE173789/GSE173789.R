
library(org.Hs.eg.db)
library(biomaRt)
library(ashr)
library(apeglm)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(ashr)


directory_path <- "E:/miRNAAA/PBMC/mRNA/12-GSE173789/GSE173789/"
MS_list <- list.files(path = directory_path, pattern = "MS", full.names = TRUE)
HC_list <- list.files(path = directory_path, pattern = "HC", full.names = TRUE)


# Read each file and combine them into one dataframe
combined_df <- data.frame()  # Initialize an empty dataframe

for (file_path in MS_list) {
  # Read the file content (adjust read function as needed)
  file_content <- read.table(file_path, header = TRUE, sep = "\t")
  
  # Append to the combined dataframe
  combined_df <-  merge(combined_df, file_content, all = TRUE)
}

for (file_path in HC_list) {
  # Read the file content (adjust read function as needed)
  file_content <- read.table(file_path, header = TRUE, sep = "\t")
  
  # Append to the combined dataframe
  combined_df <-  merge(combined_df, file_content, all = TRUE)
}

rownames(combined_df) <- combined_df$Gene_id
gen_counts <- combined_df[-1:-2] # -1 for Gene_id and -2 for NaN columns


# Create condition vector
Condition <- c(rep("MS", 22),rep("HC", 14))

# Create data frame
metadata <- data.frame(Condition )

# Assign the row names of the data frame
rownames(metadata) <- colnames(gen_counts)
all(rownames(metadata) == colnames(gen_counts))

# If "FULSE":
idx <- match(colnames(gen_counts), rownames(metadata))
ordered_metaD <- metadata[idx, ]

all(rownames(ordered_metaD) == colnames(gen_counts))


#---------------------------------------------------

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
               contrast = c("Condition", "MS","HC"),
               alpha = 0.05,lfcThreshold = 0,
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
IDs <- rownames(sig_rse)
keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)
cols <- c("SYMBOL", "GENENAME")
ensids <- IDs
gene_info <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys=ensids, columns=cols,
                                   keytype="ENSEMBL")


gene_info <- gene_info[-1,]
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
