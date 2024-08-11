library(dplyr)
library(gplots)
library(ggplot2)

processed_data <- read.table("DEsq2_results.txt", row.names = 1)
colnames(processed_data) <- c("Base.mean",	"log2FC","StdErr",
                              "Wald.Stats",	"P.value", "P.adj")

#=======================[Exporting & Visualization Of DE results]================================

#[1] Exporting

degs.res= processed_data[abs(processed_data$log2FC) >= 0  & processed_data$P.adj <= 0.05, ]
degs.res_UP=processed_data[processed_data$log2FC >= 1   & processed_data$P.adj <= 0.05,]  
degs.res_DOWN=processed_data[processed_data$log2FC  <= -1   & processed_data$P.adj <= 0.05,] 

write.table(degs.res,file = "DEGs.txt",row.names = T,col.names = T,quote = F)
write.table(rownames(degs.res_UP),file = "DEGs_UP.txt",row.names = F,col.names = F,quote = F)
write.table(rownames(degs.res_DOWN),file = "DEGs_DOWN.txt",row.names = F,col.names = F,quote = F)

#[2] Visualization
# [A] ----------- Volcano (UP/Down)-----------  

res2 <- processed_data %>%
  mutate(gene_type = case_when(log2FC >= 0  & P.adj <= 0.05 ~ "Up",
                               log2FC <= 0  & P.adj <= 0.05 ~ "Down",
                               TRUE ~ "N.s"))   
cols <- c("Up" = "#FF0000", "Down" = "#00FF00", "N.s" = "grey") 
sizes <- c("Up" = 2, "Down" = 2, "N.s" = 1) 
alphas <- c("Up" = 1, "Down" = 1, "N.s" = 0.5)

res2 %>%
  ggplot(aes(x = (log2FC),
             y = -log10(P.adj),
             fill = gene_type,    
             size = gene_type,
             alpha = gene_type)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black") + 
  geom_hline(yintercept = (0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(0),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-4, 4, 2)),       
                     limits = c(-4, 4)) 

