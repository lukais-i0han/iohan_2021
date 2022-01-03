### libraries
library('dplyr')
library('tximport')
library('DESeq2')
library('stringr')

##Metadado

Tau_individual <- read.csv('metadata/TauD35_Covariates.csv',header = T,stringsAsFactors = F)

Tau_individual <- Tau_individual %>% mutate(mutation = case_when(
  mutation == 'WT' ~ 'Control',
  TRUE ~ 'Alzheimer'
))

Tau_individual <- Tau_individual %>% mutate(age = case_when(
  age == 'Young' ~ '4_months',
  TRUE ~ '17_months'
))

Tau_individual$'Condition' <- paste(Tau_individual$mutation,Tau_individual$age,sep='_')
Tau_individual <- Tau_individual[order(Tau_individual$id),]

Tau_sex <- read.csv('metadata/tau_sex.csv',header = T, stringsAsFactors = F)
Tau_sex <- Tau_sex[order(Tau_sex$id),]

Tau_individual$'sex' <- Tau_sex$sex

Tau_individual <- Tau_individual %>% mutate(sex = case_when(
  Tau_individual$sex == '["male"]' ~ 'male',
  TRUE ~ 'female'
))

write.csv(Tau_individual,'Tau_metadado.csv')



### Import Kallisto files

files = paste0(list.files("kallisto", full.names=T), "/abundance.tsv")

tx2gene <- read.table("refs/tg2.txt", header = F, stringsAsFactors = F)

kallistoQuant <- tximport(files, type = 'kallisto',tx2gene = tx2gene[,c(1,3)])

Tau_individual <- read.csv('Tau_metadado.csv',header = T, stringsAsFactors = F, row.names = 1)
Tau_individual <- Tau_individual[order(Tau_individual$BamFileName),]

colnames(kallistoQuant[["abundance"]]) <- Tau_individual$id
colnames(kallistoQuant[["counts"]]) <- Tau_individual$id
rownames(Tau_individual) <- Tau_individual$id
### Desesq2 for gene-level

dds <- DESeqDataSetFromTximport(kallistoQuant,colData = Tau_individual,design = ~Condition)
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

vsd_TAU <- varianceStabilizingTransformation(dds,blind = F, fitType = 'local')
vsd_TAU <- assay(vsd_TAU)
saveRDS(vsd_TAU,file='results_deseq/vsd_TAU.rds')

dds <- DESeq(dds,fitType = 'local')

saveRDS(dds,file='results/Tau35.rds')

