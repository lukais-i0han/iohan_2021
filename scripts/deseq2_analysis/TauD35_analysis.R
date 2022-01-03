library('dplyr')
library('DESeq2')
library('stringr')

#### Tau deseq
Tau_deseq <- readRDS('results_deseq/Tau35.rds')

res <-results(Tau_deseq,contrast= c("Condition","Alzheimer_4_months","Control_4_months"),cooksCutoff = .99, 
              independentFiltering = T,alpha = 0.01,pAdjustMethod = "BH")

res <- res[order(res$padj),]

results <- as.data.frame(dplyr::mutate(as.data.frame(res), sig=ifelse(res$padj<0.01,
                        "FDR<0.01", "Not Sig")), row.names=rownames(res))
results<- results[order(rownames(results)),]

results <- na.omit(results)


#### Ensembl to gene symbol

genes_symbols <- read.table('refs/t2g.txt',header = F,stringsAsFactors = F, sep='\t')
genes_symbols <- genes_symbols[,c(2,3)]
colnames(genes_symbols) <- c('ENS','Gene_Symbol')
genes_symbols <- genes_symbols[!duplicated(genes_symbols$Gene_Symbol),]

row_results <- rownames(results)

genes_symbols <- genes_symbols[genes_symbols$Gene_Symbol %in% row_results,]
genes_symbols <- genes_symbols[order(genes_symbols$ENS),]

results <- results[order(rownames(results)),]

results$'Gene_Symbol' <- genes_symbols$Gene_Symbol
results$'ENSG' <- genes_symbols$ENS

HIP_17 <- results
HIP_17$'Month' <- '17'
HIP_17$'Area' <- 'Hippocampus'

###

Tau_all <- rbind(HIP_4,HIP_17)
Tau_all$'Gene_Name' <- rownames(Tau_all)

Tau_all <- Tau_all %>% mutate(DEG = case_when(
  abs(Tau_all$log2FoldChange) >= log2(1.3) & sig == 'FDR<0.01' ~ 'DEG',
  TRUE ~ 'Non_Significant'))

write.csv(Tau_all,'Tau_all.csv')

Tau_sig <- Tau_all[Tau_all$DEG != 'Non_Significant',]
Tau_sig <- Tau_sig[,c(2,6:11)]

write.csv(Tau_sig,'Tau_sig.csv')
