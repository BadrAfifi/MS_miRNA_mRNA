# Install Bioconductor Packages
BiocManager::install(c("GEOquery","GSEABase"," GOstats", "Biobase","limma"))
BiocManager::install(c("affy", "affyPLM", "frma","aroma.affymetrix"))

install.packages(c("ggplot2","gplots","pheatmap","dplyr","magrittr"))
install.packages(c("RColorBrewer", "RCurl","curl","genefilter",
                   "multtest","readr","pacman"))
BiocManager::install("hugene10stv1cdf")

library(GEOquery)
library(GSEABase)
library(GOstats)
library(Biobase)
library(limma)
library(oligo)


library(affy)
library(affyPLM)
library(frma)


library(ggplot2)
library(gplots)
library(pheatmap)
library(magrittr)# needs to be run every time you start R and want to use %>%
library(dplyr)# alternatively, this also loads %>%


library(RColorBrewer)
library(curl)
library(RCurl)
library(genefilter)
library(readr)
library(multtest)
library(pacman)

### load the CEL files

#            --- Loading the "Series" matrix data for having meta data-----

my_id <- "GSE41890"
gse <- getGEO(my_id) # Normalized data!

## check how many platforms used
length(gse)

gse1 <- gse[[1]] 


head(pData(gse1)) ## print the sample information
head(fData(gse1)) ## print the gene annotation
head(exprs(gse1)) ## print the expression data
pData(gse1)$data_processing[1] #Check the normalisation


#            --- Loading the "Affymatrix" raw data files-----

celFilesDirectory= "E:/miRNAAA/PBMC/mRNA/5-GSE41890-leukocytes/GSE41890/HC"
cels_GSE = list.celfiles(celFilesDirectory, full.names = TRUE)
cels_GSE

celFilesDirectory2= "E:/miRNAAA/PBMC/mRNA/5-GSE41890-leukocytes/GSE41890/MS"
cels_GSE2 = list.celfiles(celFilesDirectory2, full.names = TRUE)
cels_GSE2

cels_GSE3 <- c(cels_GSE,cels_GSE2)

affyData <- read.celfiles(cels_GSE3)
affyData

class(affyData)

head(exprs(affyData)) ## print the expression matrix
sampleNames(affyData)# print the sample names
head(featureNames(affyData)) ## print the probes names
head(phenoData(affyData))
annotation(affyData)

#                ----- Explotary analysis -----

#       ---[1- histogram]
cols=seq(1:length(sampleNames(affyData)))
hist(affyData ,main = "Histogram Raw data ",col=cols,target = "core")
legend(12,0.9, sampleNames(affyData),col=cols,lty=1, lwd=2,cex=0.5)

#     ---[2- box plots]
boxplot(affyData, main = "Box Plot Raw data",col=seq(1:23),target = "core")


#     ---[3- Labeling the samples]


##      ------(B) For next plots, Arranging manually "sampleInfo" which contains the type of samples
#***For your meta data, you will have to decide which columns will be useful in the analysis
# B- 
#Preparing manually according the number and location of the sample columns in the data frame

status <- c(rep("HC", 24),rep("MS", 44))
sampleInfo <- as.data.frame(status)
names(sampleInfo)[1]="Condition" # rename column 1
rownames(sampleInfo) <- (sub("\\.CEL$", "", colnames(exprs(affyData)))) 
sampleInfo



#     ---[4- Unsupervised analysis]

#***Unsupervised analysis is a good way to get an understanding of the sources of variation in the data.
#It can also identify potential outlier samples.

corMatrix <- cor(exprs(affyData),use="c") ## argument use="c" stops an error if there are any missing data points
rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix,
         annotation_col=sampleInfo, main = "Raw data"
         ,fontsize_row = 3,fontsize_col = 3)    

#     ---[5- rincipal Components Analysis (PCA)]

#**** A complementary approach is to use Principal Components Analysis (PCA).
## MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX
pca <- prcomp(t(exprs(affyData)),)
## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=Condition,label=paste("Condition", Condition))) + geom_point()

# 5.2 -Multidimensional scaling plot of distances between gene expression profiles
plotMDS(exprs(affyData), labels = sampleInfo[,"Condition"],
        gene.selection = "common", main = "Raw data")


#          ------- Data pre-processing in one step : life is so easy !--------

# threestep (background correction, normalization, summarization from probes to probesets)
# Notice: all expression measures returned by threestep are all in the log2***scale.

## S4 method for signature 'ExpressionFeatureSet'
eset <- oligo::rma(affyData,background= TRUE, 
                   normalize=TRUE, subset=NULL, target = "probeset")
#This strategy allows background subtraction, TRUE perform "RMA" background correction?
#quantile normalization and TRUE perform "quantile" normalization?
#summarization (via median-polish) #target- level of summarization (only for Exon/Gene arrays)
#subset- To be implemented.


# 1- histogram
cols=seq(1:length(sampleNames(eset)))
hist(eset ,main = "Histogram Normalized data",col=cols)
legend(12,0.9, sampleNames(eset),col=cols,lty=1, lwd=2,cex=0.5)

# 2- box plots
boxplot(eset,main = "Box Plot Normalized data ",col=seq(1:23))

# 3- Unsupervised analysis

corMatrix <- cor(exprs(eset),use="c") ## argument use="c" stops an error if there are any missing data points
rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix,
         annotation_col=sampleInfo, main = "Normalized data",fontsize_row = 3,fontsize_col = 3)    

# 4- rincipal Components Analysis (PCA)

pca <- prcomp(t(exprs(eset)),)
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=Condition,label=paste("Condition", Condition))) + geom_point()

# 4.2- Multidimensional scaling plot of distances between gene expression profiles
plotMDS(exprs(eset), labels = sampleInfo[,"Condition"],
        gene.selection = "common", main = "Normalized data")
#         ================[limma]===================

#By far the most-popular package for performing differential expression is limma.
#The design matrix is a matrix of 0 and 1s; one row for each sample and one column for each sample group.
#A 1 in a particular row and column indicates that a given sample (the row) belongs to a given group (column).

design <- model.matrix(~0+sampleInfo$Condition)
design
colnames(design) <- c("HC","MS") ## rename the columns to select them easly later 

#***Coping with outliers
#It is tempting to discard any arrays which seem to be outliers prior to differential expressions. 
#However, this is done at the expense of sample-size which could be an issue for small experiments. 
#A compromise, which has been shown to work well is to calculate **weights to define the reliability of each sample.
aw <- arrayWeights(exprs(eset),design)
aw

#***The lmFit function is used to fit the model to the data. 
#The result of which is to estimate the expression level in each of the groups that we specified.


fit <- lmFit(exprs(eset), design,weights = aw)
head(fit$coefficients)

##In order to perform the differential analysis, we have to define the contrast that we are interested in. In our case we only have two groups and one contrast of interest.
##Multiple contrasts can be defined in the makeContrasts function
# Case - Healthy (HC): gene expression of Case relative to the healthy.  
# Healthy (HC) - Case: gene expression of healthy relative to the Case.
contrasts <- makeContrasts(MS - HC, 
                           HC - MS,
                           levels=design)
## can define multiple contrasts
## e.g. makeContrasts(Group1 - Group2, Group2 - Group3,....levels=design)

fit2 <- contrasts.fit(fit, contrasts)

##Finally, apply the empirical Bayes’ step 
##to get our differential expression statistics and p-values.

fit2 <- eBayes(fit2)
head(fit2$coefficients)

topTable(fit2, coef=1) ## topTable(fit2, coef=2)
### to see the results of the second contrast (if it exists)

#If we want to know how many genes are differentially-expressed overall
#we can use the decideTests function.
table(decideTests(fit2))

## to perform adjust.method="fdr" FDR means ***Benjamini & Hochberg (1995) ("BH")***
top1 <- limma::topTable(fit2, coef=1, number=nrow(fit2), adjust.method="fdr",sort.by="none") # Cancer - Healthy
top2 <- limma::topTable(fit2, coef=2, number=nrow(fit2), adjust.method="BH",sort.by="none") #Healthy - Cancer


#   ================[Exporting & Visualization Of DE results]================

# select the Comparison items
full_results <- as.data.frame(top1) #Case-Healthy

#[1] Exporting

degs.res= full_results[abs(full_results$logFC) >= 1  & full_results$adj.P.Val <= 0.05, ]
degs.res_UP=full_results[full_results$logFC >= 1   & full_results$adj.P.Val <= 0.05,]  
degs.res_DOWN=full_results[full_results$logFC  <= -1   & full_results$adj.P.Val <= 0.05,] 

write.table(degs.res,file = "DEGs.txt",row.names = T,col.names = T,quote = F)
write.table(rownames(degs.res_UP),file = "DEGs_UP.txt",row.names = F,col.names = F,quote = F)
write.table(rownames(degs.res_DOWN),file = "DEGs_DOWN.txt",row.names = F,col.names = F,quote = F)


#[2] Visualization
# [A] ----------- Volcano (UP/Down)-----------  

res2 <- full_results %>%
  mutate(gene_type = case_when(logFC >= 2  & adj.P.Val <= 0.05 ~ "up",
                               logFC <= -1  & adj.P.Val <= 0.05 ~ "down",
                               TRUE ~ "ns"))   
cols <- c("up" = "#FF0000", "down" = "#00FF00", "ns" = "grey") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

res2 %>%
  ggplot(aes(x = (logFC),
             y = -log10(adj.P.Val),
             fill = gene_type,    
             size = gene_type,
             alpha = gene_type)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black") + 
  geom_hline(yintercept = (0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-4, 4, 2)),       
                     limits = c(-4, 4)) 

gpl <- getGEO("GPL6244")
anno <- Table(gpl)
colnames(anno)
gene_assign1<- strsplit(anno$gene_assignment, " //")
gene_symbol <- sapply(gene_assign1, function(x) x[2])

anno <- anno %>% dplyr::select ("gene_assignment")
map.df <- bind_cols(anno$ID,gene_symbol)
row.names(map.df) <- NULL
names(map.df)[1]="probe_id"
names(map.df)[2]="symbol"
