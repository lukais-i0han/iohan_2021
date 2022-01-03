library('DESeq2')
library("dplyr")
library('ggplot2')
library('ggrepel')
library("tibble")
library("stringr")

mayo_age_cer <- readRDS('results/new_results/CER/mayo_age_CER.rds')

### Function to run differential expression analysis for Human Cerebellum data from Mayo Dataset. 
### Note the we ran this analysis to Pathological Aging condition, but in the paper we decided to not included these results.

deseq_function <-function (Age,Condition) {  
  res <-results(mayo_age_cer,
                contrast = c("Diag_Age",paste0(Condition,'_',Age),paste0("Control_",Age)),
                cooksCutoff = .99, 
                independentFiltering = T,alpha = 0.01,pAdjustMethod = "BH")
  
  res <- res[order(res$padj),]
  results <- as.data.frame(dplyr::mutate(as.data.frame(res), sig=ifelse(res$padj<0.01,
                                                                        "FDR<0.01", "Not Sig")), row.names=rownames(res))
  
  results <- na.omit(results) #### 
  results$'Gene_Symbol' <- rownames(results)
  results <- results[order(results$Gene_Symbol),]
  
  results <- results %>%  mutate(DEG = case_when(
    abs(results$log2FoldChange) >= log2(1.3) & sig == 'FDR<0.01' ~ 'DEG',
    TRUE ~ 'Non_Significant'))
  
  results$'Age_Group' <- Age
  
  
  write.csv(results,file=paste0('results/',Condition,'_',Age,'_CER','.csv'))
}

Age_Group <- c('A','B','C','A','B','C','A','B')
Condition <- c('AD','AD','AD','Pathologic.Aging','Pathologic.Aging','Pathologic.Aging',
               'PSP','PSP')

for (i in 1:8) {
  deseq_function(Age_Group[i],Condition[i])
}

AD_A_CER <- read.csv('results/AD_A_CER.csv',header = T, stringsAsFactors = F,
                     row.names = 1)
AD_B_CER <- read.csv('results/AD_B_CER.csv',header = T, stringsAsFactors = F,
                     row.names = 1)
AD_C_CER <- read.csv('results/AD_C_CER.csv',header = T, stringsAsFactors = F,
                     row.names = 1)
AD_A_CER$'Group' <- 'ALZ'
AD_B_CER$'Group' <- 'ALZ'
AD_C_CER$'Group' <- 'ALZ'

PA_CER_A <- read.csv('results/Pathologic.Aging_A_CER.csv',header = T, stringsAsFactors = F,
                     row.names = 1)
PA_CER_B <- read.csv('results/Pathologic.Aging_B_CER.csv',header = T, stringsAsFactors = F,
                     row.names = 1)
PA_CER_C <- read.csv('results/Pathologic.Aging_C_CER.csv',header = T, stringsAsFactors = F,
                     row.names = 1)

PA_CER_A$'Group' <- 'PA'
PA_CER_B$'Group' <- 'PA'
PA_CER_C$'Group' <- 'PA'


PSP_A_CER <- read.csv('results/PSP_A_CER.csv',header = T, stringsAsFactors = F,
                      row.names = 1)
PSP_B_CER <- read.csv('results/PSP_B_CER.csv',header = T, stringsAsFactors = F,
                      row.names = 1)

PSP_A_CER$'Group' <- 'PSP'
PSP_B_CER$'Group' <- 'PSP'


mayo_CER_all <- rbind(AD_A_CER,AD_B_CER,AD_C_CER,PSP_A_CER,PSP_B_CER,
                      PA_CER_A,PA_CER_B,PA_CER_C)

mayo_CER_sig <- mayo_CER_all[mayo_CER_all$DEG == 'DEG',] 

write.csv(mayo_CER_all,'results/new_results/mayo_CER_all.csv') ### Table with all genes
write.csv(mayo_CER_sig,'results/new_results/mayo_CER_sig.csv') ### Table with only significant genes (FDR < 0.01)
