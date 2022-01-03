#### Import of the data

library('DESeq2')
library("dplyr")
library('ggplot2')
library('ggrepel')
library("tibble")
library("stringr")

mayo_age_TCX <- readRDS('results/mayo_age_TCX.rds')


deseq_function <-function (Age,Condition) {  
  res <-results(mayo_age_TCX,
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
  
  
  write.csv(results,file=paste0('results/',Condition,'_',Age,'_TCX','.csv'))
}

Age_Group <- c('A','B','C','A','B','C','A','B')
Condition <- c('AD','AD','AD','Pathologic.Aging','Pathologic.Aging','Pathologic.Aging',
               'PSP','PSP')

for (i in 1:8) {
  deseq_function(Age_Group[i],Condition[i])
}

AD_A_TCX <- read.csv('results/AD_A_TCX.csv',header = T, stringsAsFactors = F,
                     row.names = 1)
AD_B_TCX <- read.csv('results/AD_B_TCX.csv',header = T, stringsAsFactors = F,
                     row.names = 1)
AD_C_TCX <- read.csv('results/AD_C_TCX.csv',header = T, stringsAsFactors = F,
                     row.names = 1)
AD_A_TCX$'Group' <- 'ALZ'
AD_B_TCX$'Group' <- 'ALZ'
AD_C_TCX$'Group' <- 'ALZ'

PA_TCX_A <- read.csv('results/Pathologic.Aging_A_TCX.csv',header = T, stringsAsFactors = F,
                     row.names = 1)
PA_TCX_B <- read.csv('results/Pathologic.Aging_B_TCX.csv',header = T, stringsAsFactors = F,
                     row.names = 1)
PA_TCX_C <- read.csv('results/Pathologic.Aging_C_TCX.csv',header = T, stringsAsFactors = F,
                     row.names = 1)

PA_TCX_A$'Group' <- 'PA'
PA_TCX_B$'Group' <- 'PA'
PA_TCX_C$'Group' <- 'PA'


PSP_A_TCX <- read.csv('results/PSP_A_TCX.csv',header = T, stringsAsFactors = F,
                      row.names = 1)
PSP_B_TCX <- read.csv('results/PSP_B_TCX.csv',header = T, stringsAsFactors = F,
                      row.names = 1)

PSP_A_TCX$'Group' <- 'PSP'
PSP_B_TCX$'Group' <- 'PSP'


mayo_TCX_all <- rbind(AD_A_TCX,AD_B_TCX,AD_C_TCX,PSP_A_TCX,PSP_B_TCX,
                PA_TCX_A,PA_TCX_B,PA_TCX_C)

mayo_TCX_sig <- mayo_TCX_all[mayo_TCX_all$DEG == 'DEG',] 

write.csv(mayo_TCX_all,'results/new_results/mayo_TCX_all.csv')
write.csv(mayo_TCX_sig,'results/new_results/mayo_TCX_sig.csv')
