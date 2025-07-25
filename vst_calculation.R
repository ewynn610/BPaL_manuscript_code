#######################################
##
## Script name: VST Normalization
##
## Purpose of script: Calculate VST-normalized data with batch correction
##
## Author: Elizabeth Wynn
##
#######################################
## load the packages we will need:
library(DESeq2)

#######################################
## 01 - Get average VST expression (with batch correction)
#######################################

## Read in raw counts and meta
meta<-read.csv("data/meta.csv")
raw_counts<-read.csv("data/raw_counts.csv", row.names = 1)

## Make sure in right order
identical(meta$MARS.num, colnames(raw_counts))

## Apply VST normalization and remove batch effect
dds <- DESeqDataSetFromMatrix(countData = (raw_counts),
                              colData = meta,
                              design = ~ drug_day+batch)

vst_dat <-vst(dds, blind=F)
vst_mat <- assay(vst_dat)
mm <- model.matrix(~drug_day, colData(vst_dat))
vst_mat_batch_corrected <- limma::removeBatchEffect(vst_mat, batch=vst_dat$batch, design=mm)

## Save VST and batch corrected VST data
write.csv(vst_mat, "data/vst_dat.csv")
write.csv(vst_mat_batch_corrected, "data/batch_corrected_vst_dat.csv")

## Save DESeq2 object with batch corrected VST data
assay(vst_dat) <- vst_mat_batch_corrected

saveRDS(vst_dat, "data/vst_batch_corrected_deseq2_obj.RDS")
