### libraries
library('dplyr')
library('tximport')
library('DESeq2')
library('stringr')

## Create the metadata table to data from 5XFAD dataset

FAD_individual <- read.csv('Deseq2//metadata/UCI_5XFAD.csv',header = T,stringsAsFactors = F)
FAD_individual <- FAD_individual[!duplicated(FAD_individual$specimenID),]

FAD_individual <- FAD_individual %>% mutate(Region = case_when(
  str_detect(FAD_individual$specimenID, 'C') ~ 'CCX',
  TRUE ~ 'HIP'
)) 

FAD_individual <- FAD_individual[order(FAD_individual$individualID),]

FAD_metadata <- read.csv('Deseq2/metadata/UCI_5xFAD_individual_metadata.csv',header = T,
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

write.csv(FAD_metadata,'Deseq2/metadata/FAD_metadata_final.csv')

### Import Kallisto files

FAD_metadata <- read.csv('Deseq2/metadata/FAD_metadata_final.csv',row.names = 1,stringsAsFactors = F,
                         header = T)

files = paste0(list.files("kallisto", full.names=T), "/abundance.tsv")
files = grep(paste0(paste0("/",FAD_metadata$individualID), collapse="|"), files, value=T)

tx2gene <- read.table("Deseq2/refs/t2g.txt", header = F, stringsAsFactors = F) ### table with correspondece between transcripts and its respective genes

kallistoQuant <- tximport(files, type = 'kallisto',tx2gene = tx2gene[,c(1,3)])
colnames(kallistoQuant$abundance) <- FAD_metadata$SpecimenID
colnames(kallistoQuant$counts) <- FAD_metadata$SpecimenID

### Make table with transcript abundance to later analysis with IsoformSwitch
kallistoQuant_ISO <- tximport(files, type = 'kallisto',txOut = T)
counts_FAD5X <- data.frame(kallistoQuant_ISO$counts,stringsAsFactors = F) 
write.csv(counts_FAD5X,'counts_FAD5X.csv')

### Run DESeq2 with GMR as the design variable

dds <- DESeqDataSetFromTximport(kallistoQuant,colData = FAD_metadata,design = ~GMR)

dds <- DESeq(dds,fitType = 'local')

vsd_FAD <- varianceStabilizingTransformation(dds,blind = F, fitType = 'local')
vsd_FAD <- assay(vsd_FAD)
saveRDS(vsd_FAD,file='results_analysis/vsd_FAD.rds')

saveRDS(dds,file='Deseq2/results_analysis//5XFAD.rds')

