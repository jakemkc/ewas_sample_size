# July 18 2017
## Goal: EWAS phenotype correlation

rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
cat("\014")   # same as ctrl-L
options(max.print=3000) # default 1000, and rstudio set at 5000
options(warn=1) # default = 0. To check loop warnings



# Load -----------------------------------------------------
load("results/LIFE_df_X_Y_imputed_factor_rdy.Rdata")

# Check NA %
library(dplyr); library(magrittr)


## \\ corr ----

# select variables

lifephen <- life %>% select(SEVOLUME, SPCOUNT, TOTCNT, WHONORM, SCSADFI,SCSAHDS, PERMOT) 


colnames(lifephen) <- c("Seminal volume (ml)", "Sperm concentration (x10^6/ml)", "Total sperm count (x10^6)", 
                        "WHO morphology (%)",
                        "DNA fragmentation (%)", "High DNA stainability (%)", "Motility (%)")


library(psych)

n_phen <- corr.test(lifephen, method="spearman")$n
corrphen <- corr.test(lifephen, method="spearman")$r


halfr <- corrphen[lower.tri(corrphen, diag = FALSE)]
quantile(halfr)

## \\ heatmap ----

library(gplots)

heatmapColors <- function(numColors=16) {
    c1 <- rainbow(numColors,v=seq(0.5,1,length=numColors),s=seq(1,0.3,length=numColors),start=4/6,end=4.0001/6);
    c2 <- rainbow(numColors,v=seq(0.5,1,length=numColors),s=seq(1,0.3,length=numColors),start=1/6,end=1.0001/6);
    c3 <- c(c1,rev(c2)); 
    return(c3)
}


#


heatmap.2(corrphen,
          trace = "none", 
          margins = c(12, 12), 
          col = heatmapColors(8),  # 4 color for each +/-. It will be equally space with "breaks"' default
          srtCol = 70,  # angle of row/column labels
          srtRow = 10,  # angle of row/column labels
          cexRow = 0.8,  # row label size
          cexCol = 0.8,  # column label size
          dendrogram = "both",  # turn off row clustering. Turn off dendrogram only will still order as cluster
          Rowv = TRUE,  # row dendrogram should be reordered?
          Colv = TRUE,  # column dendrogram should be reordered?
          key.title = "",  #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          key.xlab = "",   #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          key.ylab = "",   #topleft density color key. NA to remove but resize auto; "" keeps blank wo resize
          keysize = 1,     #topleft density color key size
          key.par = list(cex = 0.45), #topleft density color text size
          symbreaks = T, # symbreaks and symkey better be speciffed together, it's the color radiation from 0 to +/- (symmetric break) 
          symkey = T)   

quartz.save("results/heatmap_EWAS_phenotype_chems.png", type = "png", device = dev.cur(), dpi = 400)
# dev.copy(png, filename="results/heatmap_r_f_chems_lipids_creat_adj.png", width=1024, height=768)
# dev.off()


