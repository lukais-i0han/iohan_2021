### libraries
library('dplyr')
library('tximport')
library('DESeq2')
library('stringr')
library('IsoformSwitchAnalyzeR')
library('readr')


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


#### Path to Kallisto files
files = paste0(list.files("kallisto", full.names=T), "/abundance.tsv")
files = grep(paste0(paste0("/",FAD_metadata$SpecimenID), collapse="|"), files, value=T)

FAD_metadata$'path' <- files

FAD_metadata[,'Month'] <- as.factor(FAD_metadata[,'Month'])
Month <- levels(FAD_metadata$Month)


#### Group of Ages

ISAR_tabs <- lapply(Month,function(M){

  myDesign <- FAD_metadata %>% filter(Month %in% M) %>% dplyr::select(SpecimenID,Group,path) %>% distinct()
  rownames(myDesign) <- rownames(myDesign) <- myDesign$SpecimenID

  kallistoQuant <- importIsoformExpression(sampleVector = myDesign$path)

  colnames(kallistoQuant$counts)[-1] <- myDesign$sampleID
  colnames(kallistoQuant$abundance)[-1] <- myDesign$sampleID
  myDesign <- myDesign[,c(1,2)]

  aSwitchList <- importRdata(
    isoformCountMatrix   = kallistoQuant$counts,
    isoformRepExpression = kallistoQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = "refs/Mus_musculus.GRCm38.99.chr_patch_hapl_scaff.gtf.gz", ### release 99
    isoformNtFasta       = "refs/Mus_musculus.GRCm38.cdna.all.fa.gz", ### release 99
    addAnnotatedORFs=T,
    showProgress = T
  )

  saveRDS(aSwitchList, paste0('results/first_analysis/Month_',M,"_aSwitchList.rds"))

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
                  addToSwitchAnalyzeRlist = T,outputPrefix = paste0('Month_',M))  

  saveRDS(aSwitchListAnalyzed,file=paste0('results/first_analysis/Month_',M,'_aSwitchListAnalyzed.rds'))

  print('This part is Done4')

})


