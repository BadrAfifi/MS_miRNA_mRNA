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

### load the CEL files

#            --- Loading the "Series" matrix data for having meta data-----

my_id <- "GSE21942"
gse <- getGEO(my_id) # Normalized data!

## check how many platforms used
length(gse)

gse1 <- gse[[1]] 


head(pData(gse1)) ## print the sample information
head(fData(gse1)) ## print the gene annotation
head(exprs(gse1)) ## print the expression data
pData(gse1)$data_processing[1] #Check the normalisation


#            --- Loading the "Affymatrix" raw data files-----

celFilesDirectory="GSE21942"
cels_GSE = list.files(celFilesDirectory, pattern = "CEL")
cels_GSE


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


##      ------For next plots, Arranging manually "sampleInfo" which contains the type of samples
#***For your meta data, you will have to decide which columns will be useful in the analysis

sampleInfo <-  pData(gse1)
colnames(sampleInfo)

## source_name_ch1 seems to contain factors we might need for the analysis. Let's pick just those columns
sampleInfo <- sampleInfo %>% dplyr::select ("disease state:ch1")
sampleInfo$`disease state:ch1` <- c(rep("HC", 15),rep("MS", 14))#Optionally, rename to shorter sample names
colnames(sampleInfo)[1] <- "Condition" #Optionally, rename to more convenient column names
sampleInfo


#     ---[4- Unsupervised analysis]

#***Unsupervised analysis is a good way to get an understanding of the sources of variation in the data.
#It can also identify potential outlier samples.

corMatrix <- cor(exprs(affyData),use="c") ## argument use="c" stops an error if there are any missing data points
rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix,
         annotation_col=sampleInfo, main = "Raw data")    

#     ---[5- rincipal Components Analysis (PCA)]

#**** A complementary approach is to use Principal Components Analysis (PCA).
## MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX
pca <- prcomp(t(exprs(affyData)),)
## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=Condition,label=paste("Condition", Condition))) + geom_point()

# 5.2
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
legend(12,0.9, sampleNames(eset),col=cols,lty=1, lwd=2,cex=0.5)

# 2- box plots
boxplot(eset,main = "Box Plot Normalized data ",col=seq(1:23))



#           ----- Aggregating the genes data

# -- Check the expression data 
# -- check values  
sum(!complete.cases(exprs(eset))) # check if there is a missing values.
sum(!is.na(exprs(eset)))# check if there is NAN 


#### mapping the probe ids into gene symbole
#### First of all, you should know what the data annotation is by annotation()
exprData1 <- bind_cols(rownames(exprs(eset)),exprs(eset))
row.names(exprData1) <- NULL
names(exprData1)[1]="probe_id" # to rename the first col to prob_id

#   ---- Getting the Gene Symbols

# B- 
anno <- fData(gse1)
colnames(anno)
anno <- anno %>% dplyr::select ("Gene Symbol")
map.df <- bind_cols(rownames(exprs(eset)),anno)
row.names(map.df) <- NULL
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
colnames(design) <- c("HC","MS") ## rename the columns to select them easly later 

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
contrasts <- makeContrasts(MS - HC,
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

#  ===========

degs.res2 <- read.table("DEGs.txt")

## Add the ' ENTREZID ' Coz it will be needed in the Pathway analyses 

mapper1 = hgu133plus2SYMBOL
mapper1
map.df1 = as.data.frame(mapper1)
head(map.df1)
symbol <- map.df[map.df1$symbol %in% rownames(degs.res2), ]

mapper2 = hgu133plus2ENTREZID
mapper2
map.df2 = as.data.frame(mapper2)
head(map.df2)
Entrez_ID <- map.df2[map.df2$probe_id %in% symbol$probe_id, ]

symbol_EntrezID <- unique(bind_cols(Entrez_ID$gene_id, symbol$symbol))
symbol_EntrezID <- symbol_EntrezID[order(symbol_EntrezID$...2), ]

# bind the IDS and symbols to DGEs file

degs.res_IDs <- bind_cols(symbol_EntrezID,degs.res)

# check if name1 and name2 columns are the same
if (identical(row.names(degs.res_IDs), degs.res_IDs$...2)) {
  print("The symbols are the same.")
} else {
  print("The columns are different.")
  different_values <- which(df$name1 != df$name2)
  print(paste("The different values are at rows:", different_values))
}

names(degs.res_IDs)[1]="IDs"
names(degs.res_IDs)[2]="symbol"
write.table(degs.res_IDs,file = "DEGs[IDs&symbol].txt",row.names = F,col.names = T,quote = F)
write.table(rownames(degs.res_IDs),file = "DEGs[names].txt",row.names = F,col.names = F,quote = F)




