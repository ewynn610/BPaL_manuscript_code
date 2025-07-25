#######################################
##
## Script name: BPaL CFU plts
##
## Purpose of script: Make CFU plots for manuscript
##
## Author: Elizabeth Wynn
##
#######################################

## load the packages we will need:
library(ggplot2)
library(scales)
library(dplyr)

#######################################
## 01 - Read in/format CFU data
#######################################

## Read in CFU
cfu=read.csv("data/bpal_cfu.csv")

## Format so that PreRx is assigned to each drug as day 0 and label # in combo
cfu=lapply(c("B", "BPaL", "P", "PaL", "BL", "BPa", "L"), function(x){
  cfu=cfu%>%filter(Group=="PreRx")%>%mutate(Group=x)
})%>%bind_rows(cfu)%>%filter(Group!="PreRx")%>%mutate(
  num_drugs=factor(case_when(Group%in% c("BL", "BPa", "PaL")~"Pairwise Combinations",
                      Group%in% c("B", "L", "P")~"Monotherapies",
                      Group =="BPaL"~"Full Regimen"),
                   levels=c("Monotherapies", 
                            "Pairwise Combinations", "Full Regimen")))

## Change factor level order for drug
cfu$Group=factor(cfu$Group, c("L", "P", "B", "PaL", "BPa", "BL", "BPaL"))

#######################################
## 02 - Get median/selected testing results
#######################################

## Calculate % change
cfu_median=cfu%>%group_by(Group, Treatment.days)%>%summarise(CFU=median(CFU, na.rm=T))

calc_perc_change=function(drug, timepoint){
  (1-(cfu_median[cfu_median$Group==drug&cfu_median$Treatment.days==timepoint,"CFU"]/
        cfu_median[cfu_median$Group==drug&cfu_median$Treatment.days==0,"CFU"]))*100
}
calc_perc_change("B", 7)
calc_perc_change("B", 21)
calc_perc_change("P", 7)
calc_perc_change("P", 21)
calc_perc_change("L", 7)
calc_perc_change("L", 21)

## Caluclate log difference
calc_log_diff=function(drug, timepoint){
  log10(cfu_median[cfu_median$Group==drug&cfu_median$Treatment.days==0,"CFU"])-
    log10(cfu_median[cfu_median$Group==drug&cfu_median$Treatment.days==timepoint,"CFU"])
}

calc_log_diff("B", 7)
calc_log_diff("B", 21)
calc_log_diff("P", 7)
calc_log_diff("P", 21)
calc_log_diff("L", 7)
calc_log_diff("L", 21)

## Get p-values for different comparisons
cfu_t.test=function(drug1, drug2, timepoint1, timepoint2){
  t.test(log10(cfu$CFU[cfu$Group==drug1&cfu$Treatment.days==timepoint1]), 
         log10(cfu$CFU[cfu$Group==drug2&cfu$Treatment.days==timepoint2]))
}

## B vs. BPa
lapply(c(2,4,7,11,14,21), function(x){
  data.frame(day=x, p.value=cfu_t.test("B","BPa", x,x)$p.value)
})%>%bind_rows()

## B vs. BPaL
lapply(c(2,4,7,11,14,21), function(x){
  data.frame(day=x, p.value=cfu_t.test("B","BPaL", x,x)$p.value)
})%>%bind_rows()
lapply(c(2,4,7,11,14,21), function(x){
  data.frame(day=x, p.value=cfu_t.test("L","L", 0,x)$p.value)
})%>%bind_rows()

##########################################
## 03 - Make plot for all drugs/regimens
##########################################

col_pal=c("L"="#ba455f","P"="#caa331", B="#6778d0", "PaL"="#c76030",
          "BPa"="#50b47b","BL"="#9750a1", "BPaL"="#919191")

ggplot(cfu, mapping=aes(as.factor(Treatment.days), CFU, group=Group, color=Group))+
  geom_point(alpha=.5, size=3)+facet_wrap(~num_drugs)+
  guides(colour = guide_legend(override.aes = list(size=6, alpha=1)))+
  stat_summary(
    fun = mean,
    geom = 'line',
    aes(group = Group, color = Group),
    size=2,show.legend = F
  )+
  theme_classic()+scale_color_manual(values=col_pal)+
  theme(legend.title=element_blank(),axis.title = element_text(size=14),strip.text = element_text(size=14),
        axis.text = element_text(size=12)
  )+
  scale_y_log10(name="Average CFU/lung",
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+xlab("Treatment Day")

#######################################
## 03 - Plot for Bpal vs. b vs. bpa
#######################################
cfu_b_bpa_bpal=cfu%>%filter(Group=="B"|Group=="BPaL"|Group=="BPa", Treatment.days!=0)

ggplot(cfu_b_bpa_bpal, aes(factor(Treatment.days), CFU, color=Group))+
  geom_boxplot(outlier.alpha = 0)+
  scale_y_log10(name="CFU/lung",
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+xlab("Treatment Day")+
  geom_point(position=position_dodge(width=.75),alpha=0.5, size=2)+theme_classic()+
  scale_color_manual(values=col_pal)+
  theme( axis.title = element_text(size=14),strip.text = element_text(size=14),
         axis.text = element_text(size=12), legend.text = element_text(size=12))+
  theme(legend.title = element_blank())
