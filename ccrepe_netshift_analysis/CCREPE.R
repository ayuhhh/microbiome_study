# ======================================================
# Title: CCREPR (data preparation for netshift analysis)
# Author: Yuhan
# Date: 07/28/2023
# ======================================================

setwd("E:\\资料积累\\代码\\肠道微生物分析\\微生物网络关联分析")

library(curatedMetagenomicData)
library(dplyr)
library(ccrepe)
library(reshape2)
library(rio)

crc_subset <- sampleMetadata %>% 
  filter(study_name=="FengQ_2015") %>%
  filter(study_condition %in% c('control', 'CRC')) %>% 
  returnSamples(dataType = "relative_abundance",
                rownames = "short")


# metadata data
meta.crc <- sampleMetadata %>%
  filter(study_name=="FengQ_2015") %>%
  filter(study_condition %in% c('control', 'CRC'))

# feature data
feat.rel.crc <- as.data.frame(assay(crc_subset))/100
feat.rel.crc <- as.data.frame(t(feat.rel.crc))

# data preparation
id.case <- meta.crc$sample_id[meta.crc$study_condition=="CRC"]
id.control <- meta.crc$sample_id[meta.crc$study_condition=="control"]

feat_case <- feat.rel.crc[id.case,]
feat_control <- feat.rel.crc[id.control,]

# case group
ccrepe.case <- ccrepe (x=feat_case, sim.score.args=list(method='spearman'))
cor.case <- ccrepe.case$sim.score
pval.case <- ccrepe.case$p.values
qval.case <- ccrepe.case$q.values

cor.case.keep <- (pval.case<0.005)*cor.case
cor.case.keep <- melt(cor.case.keep) 
colnames(cor.case.keep) <- c("feature_1", "feature_2", "cor_value")
cor.case.keep <- cor.case.keep %>% 
  subset(cor_value>0.5) 

export(cor.case.keep[, c("feature_1", "feature_2")], 
       "cor.case.keep.txt", col.names = FALSE) 

# control group
ccrepe.control <- ccrepe (x=feat_control, sim.score.args=list(method='spearman'))
cor.control <- ccrepe.control$sim.score
pval.control <- ccrepe.control$p.values
qval.control <- ccrepe.control$q.values

cor.control.keep <- (pval.control<0.005)*cor.control
cor.control.keep <- melt(cor.control.keep) 
colnames(cor.control.keep) <- c("feature_1", "feature_2", "cor_value")
cor.control.keep <- cor.control.keep %>% 
  subset(cor_value>0.5) 

export(cor.control.keep[, c("feature_1", "feature_2")], 
       "cor.control.keep.txt", col.names = FALSE) 
