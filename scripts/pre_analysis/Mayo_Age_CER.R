library('dplyr')
library('tximport')
library('DESeq2')

## Create the metadata table to data from Human cerebelum of Mayo dataset

mayo_metadado <- read.csv("refs/MayoRNAseq_RNAseq_CER_covariates.csv",
                          stringsAsFactors = F,header = T)
mayo_metadado <- mayo_metadado %>% dplyr::select(SampleID,Diagnosis,AgeAtDeath) %>% dplyr::filter(Diagnosis %in% 
                      c('Control','AD','PSP','Pathologic Aging')) %>% arrange(SampleID)
rownames(mayo_metadado) <- mayo_metadado$SampleID

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
    TRUE ~ 'Z' ## exclude patients with age below 70
  )) 

mayo_metadado$'Diag_Age' <- paste(mayo_metadado$condition,mayo_metadado$AgeAtDeath,sep='_') 


mayo_metadado <- mayo_metadado[mayo_metadado$AgeAtDeath != 'Z',]
write.csv(mayo_metadado,'results/mayo_metadado_CER.csv')

### Import of Kallisto files

files = paste0(list.files("kallisto_CER", full.names=T), "/abundance.tsv")
files = grep(paste0(paste0("/",mayo_metadado$ID), collapse="|"), files, value=T)

tx2gene <- read.table("refs/tg2.txt", header = T, stringsAsFactors = F) ### table with correspondece between transcripts and its respective genes

kallistoQuant <- tximport(files, type = 'kallisto',tx2gene = tx2gene[,c(1,3)])
colnames(kallistoQuant$abundance) <- mayo_metadado$ID
colnames(kallistoQuant$counts) <- mayo_metadado$ID

### Make table with transcript abundance to later analysis with IsoformSwitch
kallistoQuant_ISO <- tximport(files, type = 'kallisto',txOut = T)
counts_CER <- data.frame(kallistoQuant_ISO$counts,stringsAsFactors = F) 
write.csv(counts_CER,'results/new_results/counts_transcripts_CER.csv')

### Run DESeq2 with Age as the design variable

dds <- DESeqDataSetFromTximport(kallistoQuant,colData = mayo_metadado,design = ~Diag_Age)
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]
dds <- DESeq(dds)
saveRDS(dds, file= "results/new_results/mayo_age_CER.rds") 

