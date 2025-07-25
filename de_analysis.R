#######################################
##
## Script name: DE analysis
##
## Purpose of script: Run DE analysis across regimen/time conditions
##
## Author: Elizabeth Wynn
#######################################

## load up the packages we will need:
library(edgeR)
library(dplyr)

#######################################
## 01 - Read in data
#######################################

## Import raw counts
raw_counts <- read.csv("data/raw_counts.csv",row.names = 1)

## Import meta data
meta <- read.csv("data/meta.csv", stringsAsFactors = T)

## Make sure everything is in right order
identical(factor(colnames(raw_counts)), meta$MARS.num)

#######################################
## 02 - Run edgeR model
#######################################

design_mat <- model.matrix(~drug_day+batch, data=meta)
dge <- DGEList(counts = as.matrix(raw_counts))
dge<-calcNormFactors(dge)

## Run edgeR
dge <- estimateDisp(dge, design_mat)
fit_edgeR <- glmFit(dge, design_mat)

saveRDS(fit_edgeR, "output/edgeR_fit.RDS")

#######################################
## 03 - Make contrasts we're interested in
#######################################
drugs=c("B", "Pa", "L", "BL", "BPa", "PaL", "BPaL")

## Function to make contrasts
make_contrast <- function(a, b) {
  vec <- rep(0, ncol(design_mat))
  vec[which(colnames(design_mat) == a)] <- 1
  vec[which(colnames(design_mat) == b)] <- -1
  vec
}

contrast_list <- list()

# 1. Contrasts for pairs of drugs at the same timepoint
timepoints <- c(2,7,14)
for (tp in timepoints) {
  drugs_at_tp<-paste0("drug_day", drugs, "_d", tp)
  combos <- combn(drugs_at_tp, 2)
  for(i in 1:ncol(combos)) {
    name <- paste(combos[1,i], "vs", combos[2,i])
    name <- gsub("drug_day", "", name)
    contrast_list[[name]] <- make_contrast(combos[1,i], combos[2,i])
  }
}

# 2. Contrasts for  drug across timepoints
for(d in drugs) {
  drug_days <- paste0("drug_day", d, "_d", timepoints)
  combos <- combn(drug_days, 2)
  for(i in 1:ncol(combos)) {
    name <- paste(combos[1,i], "vs", combos[2,i])
    name <- gsub("drug_day", "", name)
    contrast_list[[name]] <- make_contrast(combos[1,i], combos[2,i])
  }
}
# 3. Contrasts for drugs vs PreRx_d0
pretreatment <- "PreRx_d0"
for(level in levels(meta$drug_day)) {
  if(level != pretreatment) {
    name <- paste(level, "vs", pretreatment)
    contrast_list[[name]] <- make_contrast(paste0("drug_day", level), 
                                           paste0("drug_day", pretreatment))
  }
}


## Run all contrasts
res_tabs<-lapply(contrast_list, function(x){
  res_tab<-glmLRT(fit_edgeR, contrast =x)$table
  res_tab$padj=p.adjust(res_tab$PValue, method="BH")
  res_tab
})



############## Save tables #################
saveRDS(res_tabs, "output/edgeR_res_tabs.RDS")


