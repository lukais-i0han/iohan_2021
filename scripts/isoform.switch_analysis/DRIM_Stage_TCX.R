library('stageR')
library('DRIMSeq')
library('dplyr')
library('stringr')

#### Import of ISofor Switch objects

mayo_A <- readRDS('archives/ISW/MAYO/Age_A.fullAnalysis.rds')
mayo_B <- readRDS('archives/ISW/MAYO/Age_B.fullAnalysis.rds')
mayo_C <- readRDS('archives/ISW/MAYO/Age_C.fullAnalysis.rds')

mayo_A <- mayo_A$isoformFeatures
mayo_A <- mayo_A %>% filter(condition_1 == 'AD' & condition_2 == 'Control' |
                      condition_1 == 'Control' & condition_2 == 'PSP' |
                      condition_1 == 'Control' & condition_2 == 'Pathologic_Aging',)
mayo_A <- mayo_A %>% mutate(Age_Group = case_when(
  condition_1 == 'AD' ~ 'ALZ_A',
  condition_1 == 'Control' & condition_2 == 'PSP' ~ 'PSP_A',
  TRUE ~ 'PA_A'
))


mayo_B <- mayo_B$isoformFeatures
mayo_B <- mayo_B %>% filter(condition_1 == 'AD' & condition_2 == 'Control' |
                              condition_1 == 'Control' & condition_2 == 'PSP' |
                              condition_1 == 'Control' & condition_2 == 'Pathologic_Aging')

mayo_B <- mayo_B %>% mutate(Age_Group = case_when(
  condition_1 == 'AD' ~ 'ALZ_B',
  condition_1 == 'Control' & condition_2 == 'PSP' ~ 'PSP_B',
  TRUE ~ 'PA_B'
))

mayo_C <- mayo_C$isoformFeatures
mayo_C <- mayo_C %>% filter(condition_1 == 'AD' & condition_2 == 'Control' |
                              condition_1 == 'Control' & condition_2 == 'Pathologic_Aging')
mayo_C <- mayo_C %>% mutate(Age_Group = case_when(
  condition_1 == 'AD' ~ 'ALZ_C',
  TRUE ~ 'PA_C'))

mayo_All <- rbind(mayo_A,mayo_B,mayo_C)
mayo_All <- mayo_All %>% filter(abs(dIF) > 0.05) %>% 
  dplyr::select(isoform_id,gene_id,gene_name,Age_Group)

mayo_All <- mayo_All %>% mutate(AgeAtDeath = case_when(
  str_detect(Age_Group,'_A') ~ 'A',
  str_detect(Age_Group,'_B') ~ 'B',
  TRUE ~ 'C'
))

mayo_All <- mayo_All %>% mutate(Group = case_when(
  str_detect(Age_Group,'ALZ') ~ 'ALZ',
  str_detect(Age_Group,'PSP') ~ 'PSP',
  TRUE ~ 'PA'
))

####
tx2gene <- read.table('stage_archives/refs/Human/tg2.txt',header = T)

metadado_TCX <- read.csv('stage_archives/refs/Human/new/mayo_metadado.csv',
                         stringsAsFactors = F, row.names = 1)

metadado_TCX <- metadado_TCX %>% mutate(condition = case_when(
  condition == 'AD' ~ 'ALZ',
  condition == 'Pathologic Aging' ~ 'PA',
  condition == 'Control' ~ 'Control',
  TRUE ~ 'PSP'
))

counts_all <- read.csv('stage_archives/refs/Human/new/counts_trancript_TCX.csv',
                       row.names = 1,stringsAsFactors = F)

counts_all <- counts_all[rownames(counts_all) %in% mayo_All$isoform_id,]
colnames(counts_all) <- metadado_TCX$ID

Age_Group <- c('A','B','C','A','B')
Condition <- c('PA','PA','PA','PSP','PSP')

drim_stage <- function(Age_Group,Condition){
  
    metadado <- metadado_TCX[ metadado_TCX$AgeAtDeath == Age_Group & metadado_TCX$condition == Condition 
                  | metadado_TCX$AgeAtDeath == Age_Group & metadado_TCX$condition == 'Control',]
    
    counts_local <- counts_all[,metadado$ID]
    counts_local <- counts_local[order(rownames(counts_local)),]
    
    mayo <- mayo_All[mayo_All$Group == Condition & mayo_All$AgeAtDeath == Age_Group,]

    counts_local <- counts_local[rownames(counts_local) %in% mayo$isoform_id,]
    
    tx2gene_local <- tx2gene[tx2gene$transcript %in% rownames(counts_local),]
    tx2gene_local <- tx2gene_local[order(tx2gene_local$transcript),]
    tx2gene_local <- tx2gene_local[,-2]
    colnames(tx2gene_local) <- c('feature.id','gene_id')
    
    counts_local$'feature_id' <- tx2gene_local$feature.id
    counts_local$'gene_id' <- tx2gene_local$gene_id
    
    mayo_sample <- data.frame(sample_id = metadado$ID,
                             group = metadado$condition)
    
    d <- dmDSdata(counts = counts_local,
                  samples = mayo_sample)
    
    d <- dmFilter(d, min_samps_gene_expr = length(metadado$ID), 
                  min_samps_feature_expr = 3,
                  min_gene_expr = 10, min_feature_expr = 10)
    
    design_full <- model.matrix(~condition,data=metadado)
    
    d <- dmPrecision(d,design = design_full)
    
    d <- dmFit(d, design = design_full, verbose = 1)
    
    d <- dmTest(d,coef = colnames(design_full)[2],verbose = 1)
    
    pScreen <- d@results_gene$pvalue
    names(pScreen) <- d@results_gene$gene_id
    
    pConfirmation <-matrix(d@results_feature$pvalue)
    rownames(pConfirmation) <- d@results_feature$feature_id
    
    tx2gene_local <- tx2gene_local[tx2gene_local$feature.id %in% 
                    rownames(pConfirmation),]
    
    stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
                          pScreenAdjusted = FALSE, tx2gene = tx2gene_local)
    
    stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu",
                                     alpha = 0.01,allowNA=T)
    
    padj_value <- getAdjustedPValues(stageRObj, order = F,
                                 onlySignificantGenes = F)
    
    write.csv(padj_value,
    file = paste0("stage_archives/results/new/",Condition,'_',Age_Group,".csv"))
    
    

}


for (i in 1:5) {
  drim_stage(Age_Group[i],Condition[i])
}


