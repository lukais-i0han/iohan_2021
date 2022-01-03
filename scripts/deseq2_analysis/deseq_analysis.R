### libraries
library('dplyr')
library('DESeq2')
library('stringr')


##Metadado

FAD_individual <- read.csv('metadata/UCI_5XFAD.csv',header = T,stringsAsFactors = F)
FAD_individual <- FAD_individual[!duplicated(FAD_individual$specimenID),]

FAD_individual <- FAD_individual %>% mutate(Region = case_when(
  str_detect(FAD_individual$specimenID, 'C') ~ 'CCX',
  TRUE ~ 'HIP'
)) 

FAD_individual <- FAD_individual[order(FAD_individual$individualID),]

FAD_metadata <- read.csv('metadata/UCI_5xFAD_individual_metadata.csv',header = T,
                         stringsAsFactors = F)
FAD_metadata <- FAD_metadata[FAD_metadata$individualID %in% FAD_individual$individualID,]

FAD_metadata <-FAD_metadata[rep(seq_len(nrow(FAD_metadata)), each = 2), ]

FAD_metadata <- FAD_metadata %>% mutate(Month = case_when(
  str_detect(FAD_metadata$jobNumber,'4') ~ '4',
  str_detect(FAD_metadata$jobNumber,'12') ~ '12',
  TRUE ~ '18'
  
))

FAD_metadata <- FAD_metadata %>% mutate(Group = case_when(
  str_detect(FAD_metadata$genotype,'noncarrier') ~ 'Control',
  TRUE ~ 'Alzheimer'
))


FAD_metadata <- FAD_metadata[order(FAD_metadata$individualID),]

FAD_metadata$'Region' <- FAD_individual$Region

FAD_metadata$'GMR' <- paste(FAD_metadata$Group,FAD_metadata$Month,FAD_metadata$Region,sep='_')

FAD_metadata$'SpecimenID' <- FAD_individual$specimenID
rownames(FAD_metadata) <- FAD_metadata$SpecimenID

#### write.csv(FAD_metadata,'metadata/FAD_metadata_final.csv')

#### Import of Deseq2 results

FAD_deseq <- readRDS('results_deseq/5XFAD.rds')

res <-results(FAD_deseq,contrast= c("GMR","Alzheimer_18_HIP","Control_18_HIP"),cooksCutoff = .99, 
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
genes_symbols <- genes_symbols[!duplicated(genes_symbols$ENS),]

row_results <- rownames(results)

genes_symbols <- genes_symbols[genes_symbols$ENS %in% row_results,]
genes_symbols <- genes_symbols[order(genes_symbols$ENS),]

results <- results[order(rownames(results)),]

results$'Gene_Symbol' <- genes_symbols$Gene_Symbol
results$'ENSG' <- genes_symbols$ENS

HIP_18 <- results
HIP_18$'Month' <- '18'
HIP_18$'Area' <- 'Hippocampus'

###

FAD5X_all <- rbind(CCX_4,CCX_12,CCX_18,HIP_4,HIP_12,HIP_18)

FAD5X_all <- FAD5X_all %>% mutate(DEG = case_when(
  abs(FAD5X_all$log2FoldChange) >= log2(1.3) & sig == 'FDR<0.01' ~ 'DEG',
  TRUE ~ 'Non_Significant'))

write.csv(FAD5X_all,'results_analysis/FAD5X_all.csv')

FAD5X_sig <- FAD5X_all[FAD5X_all$DEG != 'Non_Significant',]
FAD5X_sig <- FAD5X_sig[,c(2,6:12)]

write.csv(FAD5X_sig,'results_analysis/FAD5X_sig.csv')
