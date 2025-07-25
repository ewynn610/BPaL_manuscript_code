#######################################
##
## Script name: Example Clustering/Heatmap generation
##
## Purpose of script: Example of clustering/heatmap generation for average
##                    expression across timepoints. Process for B is shown
##                    A similar process was used for L, P, and B vs. BPaL
##
## Author: Elizabeth Wynn
#######################################

## load the packages we will need:
library(edgeR)
library(pheatmap)
library(ggplot2)
library(hypeR)
library(dplyr)

########################################################
## 01 - Find genes differentially expressed across time in B 
########################################################

## edgeR object in 
fit_edgeR<-readRDS("output/edgeR_fit.RDS")

## Get design matrix
design_mat<-fit_edgeR$design

## Function to make contrasts
make_contrast <- function(a, b) {
  vec <- rep(0, ncol(design_mat))
  vec[which(colnames(design_mat) == a)] <- 1
  vec[which(colnames(design_mat) == b)] <- -1
  vec
}

## Make contrasts testing for DE for B across time
contrast_names<-c(
  "B_d14_vs_PreRx_d0",
  "B_d2_vs_PreRx_d0",
  "B_d7_vs_PreRx_d0",
  "B_d2_vs_B_d14",
  "B_d7_vs_B_d14",
  "B_d7_vs_B_d2")

contrast_mat_b<-lapply(contrast_names, function(x){
  a <- paste0("drug_day",sub("_vs.*", "", x)) 
  b  <- paste0("drug_day",sub(".*vs_", "", x))
  contrast<-make_contrast(a, b)
  
  contrast
})%>%bind_cols()

## Test if any of contrasts have significant results
global_res_b=glmLRT(fit_edgeR, contrast = contrast_mat_b)$table
global_res_b$padj=p.adjust(global_res_b$PValue, method="fdr")

genes_sig_b=rownames(global_res_b)[global_res_b$padj<.05]

#######################################
## 02 - Get average average vst data (with batch correction)
#######################################

vst_batch_corrected<-read.csv("data/batch_corrected_vst_dat.csv", row.names = 1)
meta<-read.csv("data/meta.csv")

## Make sure in right order
identical(meta$MARS.num, colnames(vst_batch_corrected))

## Get average normalized data for each condition
avg_norm<-lapply(unique(meta$drug_day), function(drug_day){
  vst_sub<-vst_batch_corrected[,meta$drug_day==drug_day]
  rowMeans(vst_sub)
})%>%bind_cols()%>%data.frame()
colnames(avg_norm)<-unique(meta$drug_day)
rownames(avg_norm)<-rownames(vst_batch_corrected)

## Filter data
avg_norm_fil=avg_norm[genes_sig_b, grep("B_|PreRx",colnames(avg_norm))]%>%
  data.frame()
colnames(avg_norm_fil)[colnames(avg_norm_fil)=="PreRx_d0"] <- c("PreRx")

#######################################
## 03 - Run clustering
#######################################
## Define color scale
hmcol<-colorRampPalette(RColorBrewer::brewer.pal(10,"RdBu"))(256)

## Make starting heatmap of data
my_k=6
hm=pheatmap(as.matrix(avg_norm_fil), 
                      show_rownames = F, scale="row",
                      clustering_method="ward.D", 
                      cutree_rows = my_k, color =rev(hmcol) )

## Get original cluster assignments for each gene from cutree
## k needs to equal the number from cutree_rows
row_clusters <- cutree(hm$tree_row, k = my_k)

## Get cluster numbers from cutree to be in the order the clusters appear in the heatmap
## This doesn't happen automatically, so need to relabel
heatmap_order <- hm$tree_row$order
gene_names_in_order <- rownames(avg_norm_fil)[heatmap_order]
ordered_clusters <- row_clusters[gene_names_in_order]

## Identify the order in which clusters appear top-to-bottom
## This gives you a mapping: original cluster -> visual order
unique_clusters_in_order <- unique(ordered_clusters)

## Relabel the clusters by their visual order
## This maps original cluster numbers to 1 (top) through k (bottom)
relabel_map <- setNames(paste0("Cluster ",seq_along(unique_clusters_in_order)), unique_clusters_in_order)
relabelled_clusters <- relabel_map[as.character(row_clusters)]
## Assign gene names to new cluster assignments
names(relabelled_clusters)<-names(row_clusters)

#######################################
## 04 - Make average expression plots for each cluster
#######################################

## Loop through each cluster to make a plot for sample values over time
clust_avg_exprs_plts=lapply(paste("Cluster", 1:my_k), function(x){
  genes=names(relabelled_clusters)[relabelled_clusters==x]
  
  ## filter VST down to only genes in cluster x
  vst_fil<-vst_batch_corrected[genes, meta$MARS.num[meta$drug %in% c("B", "PreRx")]]
  
  ## Row scale the data
  vst_fil_scaled<-t(scale(t(vst_fil)))
  
  ## Get the average value across genes in cluster x for each sample
  vst_rowMeans<-data.frame(avg_exprs=colMeans(vst_fil_scaled))
  
  ## Merge with meta data
  ## Colnames of vst_rownames are MARS.nums
  vst_rowMeans=merge(vst_rowMeans, meta, by.x=0, by.y="MARS.num")
  
  ## Center around PreRx mean
  mean_PreRx=mean(vst_rowMeans$avg_exprs[ vst_rowMeans$drug=="PreRx"])
  vst_rowMeans$avg_exprs=vst_rowMeans$avg_exprs-mean_PreRx
  
  ggplot(vst_rowMeans, mapping=aes(as.factor(day), avg_exprs) )+
    geom_point(size=2.5, alpha=.55)+ 
    stat_summary(fun = mean, geom = "point", shape = "\U2013", size = 8, 
                 inherit.aes = TRUE, col = "black")+theme_classic()+
    theme(plot.title = element_text(hjust = .5, size=18))+
    ylab("Average Expression")+ggtitle(x)+xlab("Treatment Day")+
    ylim(-2.4, 2.4)+geom_hline(yintercept = 0, linetype="dashed")
  
  
})

#######################################
## 05 - Run enrichment analysis by cluster
#######################################

## Read in gene sets
cole_cats <- readRDS("data/cole_cats_list.RDS")

##  Read in curated categories
curated_sets <- readRDS("data/curated_set_list_v1.2.RDS")

## Combine curated sets and cole categories, filter out categories with <4 genes
comb_list <- c(curated_sets, cole_cats)
comb_list=comb_list[sapply(comb_list, length)>=4]

clust_enrich_res=lapply(paste("Cluster", 1:6), function(x){
  genes=names(relabelled_clusters)[relabelled_clusters==x]
  
  enrich_res <- hypeR(genes, genesets = comb_list, test="hypergeometric",
                      background=row.names(vst_batch_corrected), fdr = 1)
  
  # Making a table with relevant results, calculate percent of gene set in module
  enrich_res_tab <- enrich_res$data%>%
    arrange(fdr)%>%
    mutate(percent=(overlap/geneset)*100)%>%
    dplyr::select("Significant Genes"=overlap, "Set Size"=geneset,"Adj. P-Value"=fdr, percent)
  
})

names(clust_enrich_res)<-paste("Cluster", 1:6)
