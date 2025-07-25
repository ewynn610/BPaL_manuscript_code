#######################################
##
## Script name: Gene Set average expression analysis
##
## Purpose of script: Compare average expression of gene sets across drug conditions
##
## Author: Elizabeth Wynn
#######################################


## load up the packages we will need:
library(dplyr)
library(emmeans)
library(ggplot2)


#######################################
## 01 - Read in data
#######################################

## Read in vst data (batch corrected and not batch corrected
## Saved in vst_calculation.R
vst=read.csv("data/vst_dat.csv", row.names = 1)
vst_batch_corrected<-read.csv("data/batch_corrected_vst_dat.csv", row.names = 1)


## Read in gene sets
cole_cats <- readRDS("data/cole_cats_list.RDS")

##  Read in curated categories
curated_sets <- readRDS("data/curated_set_list_v1.2.RDS")

## Combine curated sets and cole categories
comb <- c(curated_sets, cole_cats)

## Read in meta data
meta=read.csv("data/meta.csv", stringsAsFactors = T)

#########################################################################
## 02 - Run testing of differences in avg. expression for all gene sets
#########################################################################

test_res<-lapply(names(comb), function(x){
  ## Get genes in category
  gene_names<-comb[[x]]
  gene_names<-gene_names[gene_names %in% rownames(vst)]

  ## filter VST down
  vst_fil<-vst[gene_names,]
  
  ## Row scale the data
  vst_fil_scaled<-t(scale(t(vst_fil)))
  
  ## Get the average value across genes in cluster x for each sample
  vst_rowMeans<-data.frame(avg_exprs=colMeans(vst_fil_scaled))
  
  ## Merge with meta data
  ## Colnames of vst_rownames are MARS.nums
  vst_rowMeans=merge(vst_rowMeans, meta, by.x=0, by.y="MARS.num")
  
  ## Fit models adjusting for batch
  mod<-lm(avg_exprs~drug_day + batch, data = vst_rowMeans)
  
  ## Get results for all pairwise contrast - do Tukey contrast
  res <- data.frame(emmeans(mod, pairwise ~ drug_day, adjust="none")$contrasts)
  
  ## Save category
  res$cat<-x
  res
})%>%bind_rows()

## Split results by contrast
test_res<-split(test_res, test_res$contrast)

## Do BH adjustment across p-values for gene sets in each contrast
test_res<-lapply(test_res, function(x){
  x$padj<-p.adjust(x$p.value, method = "BH")
  x
})

#######################################
## 03 - Make average expression plots
#######################################

## Selected Gene Sets
sets=c("Primary ribosomal proteins", 
       "Protein translation and modification",
       "Fatty Acid Synthases II",
       "PDIM",
       "ATP-proton motive force",
       "Cytochrome oxidase bccaa3",
       "Cytochrome oxidase bd",
       "Glycolysis",
       "Beta Oxidation",
       "TCA cycle",
       "Glyoxylate Bypass",
       "ESX1",
       "ESX3",
       "ESX4", 
       "ESX5",
       "Toxins",
       "PE and PPE families",
       "DosR",
       "IS elements Repeated sequences and Phage")

## For each gene set get avg (batch corrected) expression for each sample
avg_exprs=lapply(sets, function(set){
  ## Get genes in category
  gene_names<-comb[[set]]
  gene_names<-gene_names[gene_names %in% rownames(vst)]
  
  ## filter VST down
  vst_batch_corrected_fil<-vst_batch_corrected[gene_names,]
  
  ## Row scale the data
  vst_batch_corrected_fil_scaled<-t(scale(t(vst_batch_corrected_fil)))
  
  ## Get the average value across genes in cluster x for each sample
  vst_rowMeans<-data.frame(avg_exprs=colMeans(vst_batch_corrected_fil_scaled))
  
  ## Merge with meta data
  ## Colnames of vst_rownames are MARS.nums
  vst_rowMeans=merge(vst_rowMeans, meta, by.x=0, by.y="MARS.num")%>%
    filter(drug %in% c("PreRx", "B", "Pa", "L"))
  
  ## Get pretreatment value
  pretxt_avg_exprs=vst_rowMeans[vst_rowMeans$day==0&vst_rowMeans$drug=="PreRx",]
  
  ## Assign pretreatment as day 0 for each drug
  vst_rowMeans=lapply(c("B", "Pa", "L"), function(x){
    pretxt_avg_exprs$drug=x
    pretxt_avg_exprs
  })%>%bind_rows()%>%rbind(vst_rowMeans)%>%filter(drug!="PreRx")
  
  ## Adjust values so pretreatment mean is at 0
  mean_pretxt=mean(pretxt_avg_exprs$avg_exprs)
  vst_rowMeans$avg_exprs=vst_rowMeans$avg_exprs-mean_pretxt
  pretxt_avg_exprs$avg_exprs=pretxt_avg_exprs$avg_exprs-mean_pretxt
  vst_rowMeans$drug=factor(vst_rowMeans$drug, levels=c("L", "Pa", "B"))
  
  list(trt_exprs=vst_rowMeans, pretxt_exprs=pretxt_avg_exprs, max=max(abs(vst_rowMeans$avg_exprs)))
})
names(avg_exprs)=sets

## ylimit is the biggest absolute val from all sets
ylim=max(sapply(avg_exprs, function(x) x$max))

## Define color palette for plots
col_pal=c("L"="#ba455f","Pa"="#caa331", B="#6778d0")

## Make plot for each gene set
avg_exprs_plts=lapply(sets, function(set){
  my_avg_exprs=avg_exprs[[set]]$trt_exprs
  pretxt_avg_exprs=avg_exprs[[set]]$pretxt_exprs
  
  ggplot(my_avg_exprs, mapping=aes(as.factor(day), avg_exprs, color=drug, group=drug))+
    geom_point(size=4.5, alpha=.75)+
    stat_summary(fun = mean, geom = "line",
                 inherit.aes = TRUE, size=2, show.legend = F)+
    guides(colour = guide_legend(override.aes = list(alpha=1)))+
    theme_classic()+
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust=.5))+
    ylab("Average Scaled Expression\n(Relative to PreRx)")+#facet_wrap(~drug_new, scales="free_x")+
    xlab("Treatment Day")+
    ggtitle(set)+scale_color_manual(values=col_pal)+
    scale_x_discrete(expand=c(0,.1))+
    geom_point(data=pretxt_avg_exprs, color="#636363", size=4.5)
  
})
