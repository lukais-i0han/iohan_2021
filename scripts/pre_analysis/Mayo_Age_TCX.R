library('dplyr')
library('tximport')
library('DESeq2')
library('stringr')
library('dplyr')

##Metadado

mayo_metadado <- read.csv("refs/MayoRNAseq_RNAseq_TCX_covariates.csv",
                             stringsAsFactors = F,header = T)
mayo_metadado <- mayo_metadado %>% dplyr::select(ID,Diagnosis,AgeAtDeath) %>% dplyr::filter(Diagnosis %in% 
             c('Control','AD','PSP','Pathologic Aging')) %>% arrange(ID)
rownames(mayo_metadado) <- mayo_metadado$ID

colnames(mayo_metadado) <- c('ID','condition','AgeAtDeath')

mayo_metadado <- mayo_metadado %>%  
  mutate(AgeAtDeath = case_when(
    AgeAtDeath == '90_or_above' ~ '90',
    AgeAtDeath != '90_or_above' ~ mayo_metadado$AgeAtDeath
  )) 

mayo_metadado$AgeAtDeath <- sapply(mayo_metadado$AgeAtDeath,as.integer)

mayo_metadado <- mayo_metadado %>%  
  mutate(AgeAtDeath = case_when(
    AgeAtDeath >= 70 & AgeAtDeath <= 80 ~ 'A',
    AgeAtDeath >= 81 & AgeAtDeath <= 89 ~ 'B',
    AgeAtDeath == 90 ~ 'C',
    TRUE ~ 'Z'
  )) 

mayo_metadado$'Diag_Age' <- paste(mayo_metadado$condition,mayo_metadado$AgeAtDeath,sep='_')


mayo_metadado <- mayo_metadado[mayo_metadado$AgeAtDeath != 'Z',]
write.csv(mayo_metadado,'results/mayo_metadado.csv')
### Import Kallisto files



files = paste0(list.files("kallisto_TCX", full.names=T), "/abundance.tsv")
files = grep(paste0(paste0("/",mayo_metadado$ID), collapse="|"), files, value=T)

tx2gene <- read.table("refs/tg2.txt", header = T, stringsAsFactors = F)

kallistoQuant <- tximport(files, type = 'kallisto',tx2gene = tx2gene[,c(1,2)])
## kallistoQuant <- tximport(files, type = 'kallisto',txOut = T)
colnames(kallistoQuant$abundance) <- mayo_metadado$ID
colnames(kallistoQuant$counts) <- mayo_metadado$ID

## counts_TCX <- data.frame(kallistoQuant$counts,stringsAsFactors = F)
#write.csv(counts_TCX,'results/new_results/counts_trancript_TCX.csv')

### Desesq2 for gene-level

dds <- DESeqDataSetFromTximport(kallistoQuant,colData = mayo_metadado,design = ~Diag_Age)
keep <- rowSums(counts(dds)) >=10

dds <- dds[keep,]

dds <- DESeq(dds)
saveRDS(dds, file= "results/mayo_age_TCX.rds") 


