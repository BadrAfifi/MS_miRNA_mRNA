
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



#         --- Downloading Normalized tidy dataset as a matrix the data ----- 

my_id <- "GSE117935"
gse <- getGEO(my_id)

## check how many platforms used
length(gse)

gse2 <- gse[[1]]
gse2
class(gse2)

head(exprs(gse2))

#               ----- Explotary analysis -----

# 1- histogram
cols=seq(1:length(sampleNames(gse2)))
hist(gse2 ,main = "Histogram *Normalized data* ",col=cols)
legend(12,0.9, sampleNames(gse2),col=cols,lty=1, lwd=2,cex=0.5)

# 2- box plots
boxplot(gse2, main = "Box Plot *Normalized data*",col=seq(1:23))


##-- For next plots, *"source_name_ch1" seems to contain factors normal/case we might need for the analysis. Let's pick just those columns
# -- OR  [ Arranging Manually]
#      ------ Repreapering the dataset


ee <- as.data.frame(exprs(gse2))
dim(ee)

case_data <- ee[11:20]
cotrl_data <- ee[1:10]

exprData <- bind_cols(case_data,cotrl_data)
dim(exprData)
colnames(exprData)


##-- For next plots, Arranging manually "sampleInfo" which contains the type of samples
# -- preapering manually according the number and location of the sample columns in the data frame

samplnames <- colnames(exprData)
status <- c(rep("MS", 10),rep("HC", 10))

sampleInfo <- as.data.frame(bind_cols(samplnames,status))
names(sampleInfo)[1]="sample_names" # rename column 1
names(sampleInfo)[2]="Condition" # rename column 2
sampleInfo

# 4- Principal Components Analysis (PCA).
## MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX
pca <- prcomp(t(exprData),)
## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=Condition,label=paste("Condition", Condition))) + geom_point()




#         ================[limma]===================

#By far the most-popular package for performing differential expression is limma.
#The design matrix is a matrix of 0 and 1s; one row for each sample and one column for each sample group.
#A 1 in a particular row and column indicates that a given sample (the row) belongs to a given group (column).

design <- model.matrix(~0+sampleInfo$Condition)
design
colnames(design) <- c("Healthy","Case") ## rename the columns to select them easly later 

#***Coping with outliers
#It is tempting to discard any arrays which seem to be outliers prior to differential expressions. 
#However, this is done at the expense of sample-size which could be an issue for small experiments. 
#A compromise, which has been shown to work well is to calculate **weights to define the reliability of each sample.
aw <- arrayWeights(exprData,design)
aw

#***The lmFit function is used to fit the model to the data. 
#The result of which is to estimate the expression level in each of the groups that we specified.


fit <- lmFit(exprData, design,weights = aw)
head(fit$coefficients)

##In order to perform the differential analysis, we have to define the contrast that we are interested in. In our case we only have two groups and one contrast of interest.
##Multiple contrasts can be defined in the makeContrasts function
# Case - Healthy: gene expression of sick mice is greater than healthy mice by a value of 
# Healthy - Case: gene expression of healthy mice is greater than sick mice by a value of 
contrasts <- makeContrasts(Case - Healthy, 
                           Healthy - Case,
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
full_results <- as.data.frame(top1) #Cancer-Healthy

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




