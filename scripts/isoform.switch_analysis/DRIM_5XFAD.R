library('stageR')
library('DRIMSeq')



Age_Group <- c('CCX_4','CCX_12','CCX_18','HIP_4','HIP_12','HIP_18')

tx2gene <- read.table('stage_archives/refs/Mouse/t2g.txt',header = F)

counts_all <- readRDS('stage_archives/refs/Mouse/counts_FAD.rds')

metadado_FAD <- read.csv('stage_archives/refs/Mouse/metadado_FAD.csv',header = T, stringsAsFactors = F,
                         row.names = 1)

FAD_all <- readRDS('stage_archives/refs/Mouse/FAD_all.rds')
FAD_all$isoform_id <- str_sub(FAD_all$isoform_id,start = 1,end = 18)
FAD_all$gene_id <- str_sub(FAD_all$gene_id,start = 1,end = 18)

FAD_names <- str_sub(rownames(counts_all),start = 1,end = 18)
mart_mouse <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

results_FAD <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id",'entrezgene_id',
                                    'external_gene_name'),
                     filters = "ensembl_transcript_id", values =FAD_names,
                     mart = mart_mouse)

results_FAD <- results_FAD[!(duplicated(results_FAD$ensembl_transcript_id)),]

counts_all <- data.frame(counts_all,check.names = F)
rownames(counts_all) <- FAD_names

counts_all <- counts_all[rownames(counts_all) %in% results_FAD$ensembl_transcript_id,]

counts_all <- counts_all[order(rownames(counts_all)),]
results_FAD <- results_FAD[order(results_FAD$ensembl_transcript_id),]



## Group of Ages


stage_tabs <- lapply(Age_Group,function(Age) {
  
  
  FAD <- FAD_all[FAD_all$Group == Age,]
  
  print('Were working')
  
  metadado_FAD <- metadado_FAD[metadado_FAD$Group == Age,]
  
  counts_FAD <- counts_all[,metadado_FAD$SpecimenID]
  counts_FAD <- counts_FAD[rownames(counts_FAD) %in% FAD_all$isoform_id,]
  
  FAD_id <- results_FAD[results_FAD$ensembl_transcript_id %in% rownames(counts_FAD),]
  FAD_id <- FAD_id[,c(2,1)]
  colnames(FAD_id) <- c('feature_id','gene_id')
  
  counts_FAD$'feature_id' <- FAD_id$'feature_id'
  counts_FAD$'gene_id' <- FAD_id$'gene_id'
  
  FAD_samples <- data.frame(sample_id = metadado_FAD$SpecimenID,
                              group = metadado_FAD$GMR)
  
  d <- dmDSdata(counts = counts_FAD,samples = FAD_samples)
  
  d <- dmFilter(d, min_samps_gene_expr = length(metadado_FAD$individualID), min_samps_feature_expr = 3,
                min_gene_expr = 10, min_feature_expr = 10)
  design_full <- model.matrix(~ GMR,data=metadado_FAD)
  
  d <- dmPrecision(d,design = design_full)
  
  d <- dmFit(d, design = design_full, verbose = 1)
  
  d <- dmTest(d,coef = colnames(design_full)[2],verbose = 1)
  
  pScreen <- d@results_gene$pvalue
  names(pScreen) <- d@results_gene$gene_id
  
  pConfirmation <-matrix(d@results_feature$pvalue)
  rownames(pConfirmation) <- d@results_feature$feature_id
  
  stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
                        pScreenAdjusted = FALSE, tx2gene = FAD_id)
  
  stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu",
                                   alpha = 0.01,allowNA=T)
  getSignificantGenes(stageRObj)
  getSignificantTx(stageRObj)
  
  padj <- getAdjustedPValues(stageRObj, order = TRUE,
                               onlySignificantGenes = T)
  
  write.csv(padj,file = paste0("stage_archives/results/Mouse/stage_FAD_",Age,".csv"))
  
  print("Over")
})

