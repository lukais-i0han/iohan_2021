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


### Import Kallisto files

files = paste0(list.files("kallisto", full.names=T), "/abundance.h5")

tx2gene <- read.table("refs/tg2.txt", header = F, stringsAsFactors = F)

Tau_individual$'path' <- files

Tau_individual[,'age'] <- as.factor(Tau_individual[,'age'])
age <- levels(Tau_individual$age)

#### Group of Ages

ISAR_tabs <- lapply(age,function(A){
  
  myDesign <- Tau_individual %>% filter(age %in% M) %>% dplyr::select(id,Condition,path) %>% distinct()
  rownames(myDesign) <- rownames(myDesign) <- myDesign$id
  
  kallistoQuant <- importIsoformExpression(sampleVector = myDesign$path)
  
  colnames(kallistoQuant$counts)[-1] <- myDesign$sampleID
  colnames(kallistoQuant$abundance)[-1] <- myDesign$sampleID
  myDesign <- myDesign[,c(1,2)]
  colnames(myDesign) <- c('sampleID','condition')
  
  aSwitchList <- importRdata(
    isoformCountMatrix   = kallistoQuant$counts,
    isoformRepExpression = kallistoQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = "refs/Mus_musculus.GRCm38.99.chr_patch_hapl_scaff.gtf.gz", ### release 99
    isoformNtFasta       = "refs/Mus_musculus.GRCm38.cdna.all.fa.gz", ### release 99
    addAnnotatedORFs=T,
    showProgress = T
  )
  
  saveRDS(aSwitchList, paste0('results/first_analysis/age_',A,"_aSwitchList.rds"))
  
  print('This part is Done1')
  
  aSwitchListFiltered <- preFilter(
    switchAnalyzeRlist = aSwitchList,
    geneExpressionCutoff = 10,
    isoformExpressionCutoff = 3,
    removeSingleIsoformGenes = T
  )
  
  print('This part is Done2')
  
  aSwitchListAnalyzed <- isoformSwitchTestDRIMSeq(
    switchAnalyzeRlist = aSwitchListFiltered,
    reduceFurtherToGenesWithConsequencePotential=FALSE,
    dIFcutoff=0.05,
    reduceToSwitchingGenes=TRUE
  )
  
  aSwitchListAnalyzed <- analyzeORF(
    aSwitchListAnalyzed,
    orfMethod = "longest",
    showProgress= T
  )
  
  
  aSwitchListAnalyzed <- extractSequence(aSwitchListAnalyzed,writeToFile=F,quiet=T)
  print('This part is Done3')
  
  #Extract Sequence
  
  extractSequence(aSwitchListAnalyzed,pathToOutput = 'results/first_analysis/fasta/',
                  addToSwitchAnalyzeRlist = T,outputPrefix = paste0('age_',A))  
  
  saveRDS(aSwitchListAnalyzed,file=paste0('results/first_analysis/age_',A,'_aSwitchListAnalyzed.rds'))
  
  print('This part is Done4')
  
})


