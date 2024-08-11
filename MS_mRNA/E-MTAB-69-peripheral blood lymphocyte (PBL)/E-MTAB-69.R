#                --- Installing Packages and Loading -----

install.packages("BiocManager")

# Install Bioconductor Packages
BiocManager::install(c("GEOquery","GSEABase"," GOstats", "Biobase","limma"))
BiocManager::install(c("affy", "affyPLM", "frma"))

install.packages(c("ggplot2","gplots","pheatmap","dplyr","magrittr"))
install.packages(c("RColorBrewer", "RCurl","curl","genefilter",
                   "multtest","readr","pacman"))


library(GEOquery)
library(GSEABase)
library(GOstats)
library(Biobase)
library(limma)

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


#            --- Loading the "Affymatrix" raw data files-----

celFilesDirectoryMS= "E:/miRNAAA/PBMC/mRNA/8-E-MTAB-69-peripheral blood lymphocyte (PBL)/MS"
cels_GSEMS = list.celfiles(celFilesDirectoryMS, full.names = TRUE)
cels_GSEMS

celFilesDirectoryND= "E:/miRNAAA/PBMC/mRNA/8-E-MTAB-69-peripheral blood lymphocyte (PBL)/ND"
cels_GSEND = list.celfiles(celFilesDirectoryND, full.names = TRUE)
cels_GSEND

celFilesDirectory <- c(celFilesDirectoryMS,celFilesDirectoryND)

affyData <- ReadAffy(celfile.path=celFilesDirectory)
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
hist(affyData ,main = "Histogram Raw data ",col=cols)
legend(12,0.9, sampleNames(affyData),col=cols,lty=1, lwd=2,cex=0.5)

#     ---[2- box plots]
boxplot(affyData, main = "Box Plot Raw data",col=seq(1:23))


#     ---[3- Labeling the samples]

# B- 
#Preparing manually according the number and location of the sample columns in the data frame

status <- c(rep("MS", 24),rep("ND", 18))
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
         annotation_col=sampleInfo, main = "Raw data",fontsize_row = 3,fontsize_col = 3)    

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

eset = threestep(affyData,
                 background.method = "IdealMM",
                 normalize.method = "quantile",
                 summary.method = "median.polish")
# 1- histogram
cols=seq(1:length(sampleNames(eset)))
hist(eset,main = "Histogram Normalized data",col=cols)
legend(12,10, sampleNames(eset),col=cols,lty=1, lwd=2,cex=0.5)

# 2- box plots
boxplot(eset,main = "Box Plot Normalized data ",col=seq(1:23),fontsize_col = 1)

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


#           ----- Aggregating the genes data

# -- Check the expression data 
# -- check values  
sum(!complete.cases(exprs(eset))) # check if there is a missing values.
sum(is.na(exprs(eset)))# check if there is NAN 


#### mapping the probe ids into gene symbole
#### First of all, you should know what the data annotation is by annotation()
exprData1 <- bind_cols(rownames(exprs(eset)),exprs(eset))
row.names(exprData1) <- NULL
names(exprData1)[1]="probe_id" # to rename the first col to prob_id

#   ---- Getting the Gene Symbols
# A-
annotation(affyData)
# BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
hgu133plus2()
mapper = hgu133plus2SYMBOL
mapper
map.df = as.data.frame(mapper)
head(map.df) # two columns together "probe_id" and "symbol"
names(map.df)[1]="probe_id"
names(map.df)[2]="symbol"


#    ---- Then:
# merge the two data frames to have the symbole annotation in the data object
exprData2 <- merge(exprData1,map.df,by="probe_id",all.x=T)

# do i need the probe id again?  no, then let's drop it
exprData2=exprData2[,-1]

# remove nulls : some probes were not mapped to any gene symbol
sum(is.na(exprData2$symbol))
exprData2 <- exprData2[ ! is.na(exprData2$symbol),]
dim(exprData2)

# check duplciation of of gene symbols?  
sum(duplicated(exprData2$symbol))


### yes .. why ? probesets?  solutions : aggregation
agg_exprData= aggregate(exprData2, list(exprData2$symbol),FUN=mean)
dim(agg_exprData)

rownames(agg_exprData) <- agg_exprData[,1]# to rename the gene symbols as rownames
agg_exprData <- agg_exprData[,-1] # remove the gene symbols as a column
agg_exprData$symbol <- NULL # remove the gene symbols as a column


#         ================[limma]===================

#By far the most-popular package for performing differential expression is limma.
#The design matrix is a matrix of 0 and 1s; one row for each sample and one column for each sample group.
#A 1 in a particular row and column indicates that a given sample (the row) belongs to a given group (column).

design <- model.matrix(~0+sampleInfo$Condition)
design
colnames(design) <- c("MS","ND") ## rename the columns to select them easly later 

#***Coping with outliers
#It is tempting to discard any arrays which seem to be outliers prior to differential expressions. 
#However, this is done at the expense of sample-size which could be an issue for small experiments. 
#A compromise, which has been shown to work well is to calculate **weights to define the reliability of each sample.
aw <- arrayWeights(agg_exprData,design)
aw

#***The lmFit function is used to fit the model to the data. 
#The result of which is to estimate the expression level in each of the groups that we specified.


fit <- lmFit(agg_exprData, design,weights = aw)
head(fit$coefficients)

##In order to perform the differential analysis, we have to define the contrast that we are interested in. In our case we only have two groups and one contrast of interest.
##Multiple contrasts can be defined in the makeContrasts function
# Case - Healthy (HC): gene expression of Case relative to the healthy.  
# Healthy (HC) - Case: gene expression of healthy relative to the Case.
contrasts <- makeContrasts(MS - ND,
                           levels=design)

## can define multiple contrasts
## e.g. makeContrasts(Group1 - Group2, Group2 - Group3,....levels=design)

fit2 <- contrasts.fit(fit, contrasts)

##Finally, apply the empirical Bayes’ step 
##to get our differential expression statistics and p-values.
fit2 <- eBayes(fit2)
head(fit2$coefficients)

#If we want to know how many genes are differentially-expressed overall
#we can use the decideTests function.
table(decideTests(fit2))

topTable(fit2, coef=1) ## topTable(fit2, coef=2) ### to see the results of the second contrast (if it exists)


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
  mutate(gene_type = case_when(logFC >= 1  & adj.P.Val <= 0.05 ~ "up",
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


#[2] Visualization
# [B] ----------- HeatMap(DEGs)-----------  

# [ Selecting the set]
full_results <- topTable(fit2, coef=1, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"ID")
# [1]
my_genes <- rownames(degs.res)
# [2] # take care of Gene.Symbol maybe Symbol!!! depends on Annotation
ids_of_interest <-  filter(full_results,Gene.Symbol %in% my_genes) %>% 
  pull(ID)
# [3] # take care of Gene.Symbol maybe Symbol!!! depends on Annotation
gene_names <-  filter(full_results, ID  %in% my_genes) %>% 
  pull(ID)
# [4]
ex <- agg_exprData
gene_matrix <- ex[ids_of_interest,]
# [5] from 'sampleInfo'

annotation <- sampleInfo
rownames(annotation)= colnames(ex)
names(annotation)[1] <- "Sample"

# Specify colors
Sample = c("lightgreen", "navy")
names(Sample) = c("HC", "MS")
ann_colors = list(Sample = Sample)

# Convort the color use sampleInfo instead of annotation in the following code

pheatmap(gene_matrix,labels_row = ids_of_interest, annotation = sampleInfo, annotation_colors = ann_colors, scale="row",fontsize_row = 4,fontsize_col = 5)
pheatmap(gene_matrix,labels_row = gene_names, annotation = annotation, annotation_colors = ann_colors, scale="row",fontsize_row = 4,fontsize_col = 5)

heatmap.2(gene_matrix,labels_row = gene_names, cexCol=0.7,col = rev(redblue(26)), scale = "row")

sessionInfo()