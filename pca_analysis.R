#######################################
##
## Script name: Run PCA Analysis
##
## Purpose of script: Run PCA and make analysis plots
##
## Author: Elizabeth Wynn
##
#######################################

## load the packages we will need:
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(ggConvexHull)
library(dplyr)

#######################################
## 01 - Read in data
#######################################
vst_batch_corrected<-readRDS("data/vst_batch_corrected_deseq2_obj.RDS")

#######################################
## 02 - Run PCA
#######################################
pca_dat_bydrug=plotPCA(vst_batch_corrected, intgroup=c("drug_day", "drug", "day"),
                       returnData=T)

percentVar <- round(100 * attr(pca_dat_bydrug, "percentVar")) 

#######################################
## 03 - Make PCA plots
#######################################

## Color palette
col_pal=c("B"="#6778d0","Pa"="#caa331", L="#ba455f", "PaL"="#c76030",
          "BPa"="#50b47b","BL"="#9750a1", "PreRx" = "black", "BPaL"="#919191")

## Format data
pca_dat_bydrug$drug<-factor(pca_dat_bydrug$drug, 
                            levels=c("PreRx","B", "Pa", "L", 
                                     "BPa", "BL", "PaL", "BPaL"))

pca_dat_bydrug$day=factor(paste0("Day ",pca_dat_bydrug$day),
                          levels=c("Day 0", "Day 2", "Day 7", "Day 14"))


ggplot(pca_dat_bydrug, aes(x = PC1, y = PC2, color = drug)) +
  geom_point(size =3.5, mapping=aes(shape=factor(day))) + 
  scale_shape_manual(values=15:18)+
  scale_color_manual("",values = col_pal)+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic()+
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size=18),
        axis.text = element_text(size=15),
        legend.title = element_blank())

#######################################
## 04 - Make plots of individual genes
#######################################

## Filter down to monotherapies and preRx
pca_dat_monotherapies=pca_dat_bydrug%>%
  filter(drug %in% c("B", "Pa", "L"))

pca_dat_preRx=pca_dat_bydrug%>%
  filter(drug %in% c("PreRx"))

## Assign PreRx to each drug as day 0
PreRx_dat_drug=lapply(c("B","L","Pa"), function(x){
  pca_dat_preRx$drug=x
  pca_dat_preRx
})%>%bind_rows()

pca_dat_monotherapies=bind_rows(pca_dat_monotherapies, PreRx_dat_drug)%>%
  mutate(drug=factor(drug, levels=c("L", "Pa", "B")))

## Mean PCA for labels
mean_pc=pca_dat_monotherapies%>%
  group_by(drug_day, drug, day)%>%
  summarise(mean_pc1=mean(PC1), mean_pc2=mean(PC2))%>%arrange(day)

mean_pc_and_pts=data.frame(drug_day=NA, drug=pca_dat_monotherapies$drug, day="", 
                           mean_pc1=pca_dat_monotherapies$PC1,
                           mean_pc2=pca_dat_monotherapies$PC2)%>%rbind(mean_pc)

ggplot(pca_dat_monotherapies, aes(x = PC1, y = PC2, color = drug)) + 
  geom_convexhull(aes(group=day, fill=drug),color="white",alpha=.5,show.legend = F)+
  geom_point(size =3.5) + 
  scale_color_manual("",values = col_pal)+
  scale_fill_manual("",values = col_pal)+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  facet_wrap(~drug, ncol = 1)+
  theme_classic()+
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size=24),
        axis.text = element_text(size=18),
        legend.position = "none")+
  ggrepel::geom_label_repel(data=mean_pc_and_pts, aes(x=mean_pc1, y=mean_pc2,
                                                      label=day),
                            label.padding=unit(0.15, "lines"),size=4.5, min.segment.length = 0,
                            seed=24, force = 30, color="black")

#######################################################
## 05 - Make plots for pairwise comparisons at day 14
#######################################################

## Make data frame of subsets to make plots
drug_combs=data.frame(drug1=c("L", "Pa", "L"), drug2=c("B", "B", "Pa"))%>%
  mutate(combo=paste0(drug2, drug1))

## Make plots with pairwise combinations and their individual components
pairwise_comb_plots=apply(drug_combs, 1, function(x){
  
  pca_dat_fil=pca_dat_bydrug%>%filter(drug %in% c(x, "PreRx"), day %in% c("Day 0", "Day 14"))%>%
    mutate(drug=factor(drug, levels=c(x,"PreRx"), labels = c(x, "PreRx")))
  
  ggplot(pca_dat_fil, aes(x = PC1, y = PC2, color = drug)) + 
    geom_convexhull(aes(group=drug_day, fill=drug),alpha=.5,show.legend = F)+
    geom_point(size =3.5) + 
    scale_color_manual("",values = col_pal)+
    scale_fill_manual("",values = col_pal)+
    xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
    ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
    theme_classic()+
    theme(legend.title = element_blank())+
    xlim(-34,42)+ylim(-22,15)
  
})

#######################################
## 06 - Make plot subset with B and BPaL
#######################################

pca_dat_fil=pca_dat_bydrug%>%filter(drug %in% c("B","BPaL", "PreRx"), day %in% c("Day 0", "Day 14"))%>%
  mutate(drug=factor(drug, levels=c("B","BPaL","PreRx"), labels = c("B","BPaL", "PreRx")))

ggplot(pca_dat_fil, aes(x = PC1, y = PC2, color = drug)) + 
  geom_convexhull(aes(group=drug_day, fill=drug),alpha=.5,show.legend = F)+
  geom_point(size =3.5) + 
  scale_color_manual("",values = col_pal)+
  scale_fill_manual("",values = col_pal)+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  theme_classic()+
  theme(legend.title = element_blank())
