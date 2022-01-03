library('gprofiler2')


F_DEG_DTU <- readRDS('NEW_GPROFILER/refs/F_DEG_DTU.rds')



### DEG

F_DEG <- F_DEG_DTU[F_DEG_DTU$Class == 'DEG',]

DEG_FAD <- gost(query = list(
  'CCX[4]' = F_DEG[F_DEG$Area == 'Cortex' & F_DEG$Age_Group == '[4]'|
             F_DEG$Area == 'Cortex' & F_DEG$Age_Group == '[4-12]'|
             F_DEG$Area == 'Cortex' & F_DEG$Age_Group == '[4-18]',1],
  'CCX[12]' = F_DEG[F_DEG$Area == 'Cortex' & F_DEG$Age_Group == '[12]'|
              F_DEG$Area == 'Cortex' & F_DEG$Age_Group == '[4-12]'|
              F_DEG$Area == 'Cortex' & F_DEG$Age_Group == '[12-18]',1],
  'CCX[18]' = F_DEG[F_DEG$Area == 'Cortex' & F_DEG$Age_Group == '[18]'|
              F_DEG$Area == 'Cortex' & F_DEG$Age_Group == '[4-18]' |
              F_DEG$Area == 'Cortex' & F_DEG$Age_Group == '[12-18]',1],
  'CCX[4.12.18]' = F_DEG[F_DEG$Area == 'Cortex' & F_DEG$Age_Group == '[4-12-18]'|
                   F_DEG$Area == 'Cortex' & F_DEG$Age_Group == '[4-12]' |
                   F_DEG$Area == 'Cortex' & F_DEG$Age_Group == '[4-18]' |
                   F_DEG$Area == 'Cortex' & F_DEG$Age_Group == '[12-18]',1],
  'HIP[4]' = F_DEG[F_DEG$Area == 'Hippocampus' & F_DEG$Age_Group == '[4]'|
             F_DEG$Area == 'Hippocampus' & F_DEG$Age_Group == '[4-12]' |
             F_DEG$Area == 'Hippocampus' & F_DEG$Age_Group == '[4-18]',1],
  'HIP[12]' = F_DEG[F_DEG$Area == 'Hippocampus' & F_DEG$Age_Group == '[12]'|
              F_DEG$Area == 'Hippocampus' & F_DEG$Age_Group == '[4-12]'|
              F_DEG$Area == 'Hippocampus' & F_DEG$Age_Group == '[12-18]',1],
  'HIP[18]' = F_DEG[F_DEG$Area == 'Hippocampus' & F_DEG$Age_Group == '[18]'|
              F_DEG$Area == 'Hippocampus' & F_DEG$Age_Group == '[4-18]' |
              F_DEG$Area == 'Hippocampus' & F_DEG$Age_Group == '[12-18]',1],
  'HIP[4.12.18]' = F_DEG[F_DEG$Area == 'Hippocampus' & F_DEG$Age_Group == '[4-12-18]'|
                   F_DEG$Area == 'Hippocampus' & F_DEG$Age_Group == '[4-12]' |
                   F_DEG$Area == 'Hippocampus' & F_DEG$Age_Group == '[4-18]' |
                   F_DEG$Area == 'Hippocampus' & F_DEG$Age_Group == '[12-18]',1]),
  organism = 'mmusculus',evcodes = T,
  significant = T,
  multi_query = F,sources = c('GO:BP','GO:MF','GO:CC'),
  user_threshold = 0.01,
  correction_method = 'fdr')



DEG_FAD <- DEG_FAD$result[,c("query", "term_id", "term_name", "p_value", "intersection",
                             'precision','intersection_size')]
colnames(DEG_FAD) <-  c("query", "GO.ID", "Description", "p.Val",'Genes','precision',
                        'intersection_size')
DEG_FAD$FDR <- DEG_FAD$p.Val
DEG_FAD$"Class" <- 'DEG'
DEG_FAD$'Phenotype' <- "+1"
DEG_FAD <- DEG_FAD[DEG_FAD$precision >0.03,]



DEG_FAD %>% group_by(query) %>%
  group_walk(~
               write.table(data.frame(.x[,c("GO.ID", "Description", "p.Val", "FDR",'Phenotype',"Genes")]), 
                           file = paste0("NEW_GPROFILER/FAD/DEG/", unique(.y$query), "_gem.txt"),
                           sep = "\t", quote = F, row.names = F))


### DTU


F_DTU <- F_DEG_DTU[F_DEG_DTU$Class == 'DTU',]

DTU_FAD <- gost(query = list(
  'CCX[4]' = F_DTU[F_DTU$Area == 'Cortex' & F_DTU$Age_Group == '[4]'|
                     F_DTU$Area == 'Cortex' & F_DTU$Age_Group == '[4-12]'|
                     F_DTU$Area == 'Cortex' & F_DTU$Age_Group == '[4-18]',1],
  'CCX[12]' = F_DTU[F_DTU$Area == 'Cortex' & F_DTU$Age_Group == '[12]'|
                      F_DTU$Area == 'Cortex' & F_DTU$Age_Group == '[4-12]'|
                      F_DTU$Area == 'Cortex' & F_DTU$Age_Group == '[12-18]',1],
  'CCX[18]' = F_DTU[F_DTU$Area == 'Cortex' & F_DTU$Age_Group == '[18]'|
                      F_DTU$Area == 'Cortex' & F_DTU$Age_Group == '[4-18]' |
                      F_DTU$Area == 'Cortex' & F_DTU$Age_Group == '[12-18]',1],
  'CCX[4.12.18]' = F_DTU[F_DTU$Area == 'Cortex' & F_DTU$Age_Group == '[4-12-18]'|
                           F_DTU$Area == 'Cortex' & F_DTU$Age_Group == '[4-12]' |
                           F_DTU$Area == 'Cortex' & F_DTU$Age_Group == '[4-18]' |
                           F_DTU$Area == 'Cortex' & F_DTU$Age_Group == '[12-18]',1],
  'HIP[4]' = F_DTU[F_DTU$Area == 'Hippocampus' & F_DTU$Age_Group == '[4]'|
                     F_DTU$Area == 'Hippocampus' & F_DTU$Age_Group == '[4-12]' |
                     F_DTU$Area == 'Hippocampus' & F_DTU$Age_Group == '[4-18]',1],
  'HIP[12]' = F_DTU[F_DTU$Area == 'Hippocampus' & F_DTU$Age_Group == '[12]'|
                      F_DTU$Area == 'Hippocampus' & F_DTU$Age_Group == '[4-12]'|
                      F_DTU$Area == 'Hippocampus' & F_DTU$Age_Group == '[12-18]',1],
  'HIP[18]' = F_DTU[F_DTU$Area == 'Hippocampus' & F_DTU$Age_Group == '[18]'|
                      F_DTU$Area == 'Hippocampus' & F_DTU$Age_Group == '[4-18]' |
                      F_DTU$Area == 'Hippocampus' & F_DTU$Age_Group == '[12-18]',1],
  'HIP[4.12.18]' = F_DTU[F_DTU$Area == 'Hippocampus' & F_DTU$Age_Group == '[4-12-18]'|
                           F_DTU$Area == 'Hippocampus' & F_DTU$Age_Group == '[4-12]' |
                           F_DTU$Area == 'Hippocampus' & F_DTU$Age_Group == '[4-18]' |
                           F_DTU$Area == 'Hippocampus' & F_DTU$Age_Group == '[12-18]',1]),
  organism = 'mmusculus',evcodes = T,
  significant = T,
  multi_query = F,sources = c('GO:BP','GO:MF','GO:CC'),
  user_threshold = 0.01,
  correction_method = 'fdr')



DTU_FAD <- DTU_FAD$result[,c("query", "term_id", "term_name", "p_value", "intersection",
                             'precision','intersection_size')]
colnames(DTU_FAD) <-  c("query", "GO.ID", "Description", "p.Val",'Genes','precision',
                        'intersection_size')
DTU_FAD$FDR <- DTU_FAD$p.Val
DTU_FAD$"Class" <- 'DEG'
DTU_FAD$'Phenotype' <- "+1"
DTU_FAD <- DTU_FAD[DTU_FAD$precision >0.03,]

DTU_FAD %>% group_by(query) %>%
  group_walk(~
               write.table(data.frame(.x[,c("GO.ID", "Description", "p.Val", "FDR",'Phenotype',"Genes")]), 
                           file = paste0("NEW_GPROFILER/FAD/DTU/", unique(.y$query), "_gem.txt"),
                           sep = "\t", quote = F, row.names = F))

### DEG-DTU


DEG_DTU_FAD <- gost(query = list(
  'CCX[4]' = F_DEG_DTU[F_DEG_DTU$Area == 'Cortex' & F_DEG_DTU$Age_Group == '[4]'|
                     F_DEG_DTU$Area == 'Cortex' & F_DEG_DTU$Age_Group == '[4-12]'|
                     F_DEG_DTU$Area == 'Cortex' & F_DEG_DTU$Age_Group == '[4-18]',1],
  'CCX[12]' = F_DEG_DTU[F_DEG_DTU$Area == 'Cortex' & F_DEG_DTU$Age_Group == '[12]'|
                      F_DEG_DTU$Area == 'Cortex' & F_DEG_DTU$Age_Group == '[4-12]'|
                      F_DEG_DTU$Area == 'Cortex' & F_DEG_DTU$Age_Group == '[12-18]',1],
  'CCX[18]' = F_DEG_DTU[F_DEG_DTU$Area == 'Cortex' & F_DEG_DTU$Age_Group == '[18]'|
                      F_DEG_DTU$Area == 'Cortex' & F_DEG_DTU$Age_Group == '[4-18]' |
                      F_DEG_DTU$Area == 'Cortex' & F_DEG_DTU$Age_Group == '[12-18]',1],
  'CCX[4.12.18]' = F_DEG_DTU[F_DEG_DTU$Area == 'Cortex' & F_DEG_DTU$Age_Group == '[4-12-18]'|
                           F_DEG_DTU$Area == 'Cortex' & F_DEG_DTU$Age_Group == '[4-12]' |
                           F_DEG_DTU$Area == 'Cortex' & F_DEG_DTU$Age_Group == '[4-18]' |
                           F_DEG_DTU$Area == 'Cortex' & F_DEG_DTU$Age_Group == '[12-18]',1],
  'HIP[4]' = F_DEG_DTU[F_DEG_DTU$Area == 'Hippocampus' & F_DEG_DTU$Age_Group == '[4]'|
                     F_DEG_DTU$Area == 'Hippocampus' & F_DEG_DTU$Age_Group == '[4-12]' |
                     F_DEG_DTU$Area == 'Hippocampus' & F_DEG_DTU$Age_Group == '[4-18]',1],
  'HIP[12]' = F_DEG_DTU[F_DEG_DTU$Area == 'Hippocampus' & F_DEG_DTU$Age_Group == '[12]'|
                      F_DEG_DTU$Area == 'Hippocampus' & F_DEG_DTU$Age_Group == '[4-12]'|
                      F_DEG_DTU$Area == 'Hippocampus' & F_DEG_DTU$Age_Group == '[12-18]',1],
  'HIP[18]' = F_DEG_DTU[F_DEG_DTU$Area == 'Hippocampus' & F_DEG_DTU$Age_Group == '[18]'|
                      F_DEG_DTU$Area == 'Hippocampus' & F_DEG_DTU$Age_Group == '[4-18]' |
                      F_DEG_DTU$Area == 'Hippocampus' & F_DEG_DTU$Age_Group == '[12-18]',1],
  'HIP[4.12.18]' = F_DEG_DTU[F_DEG_DTU$Area == 'Hippocampus' & F_DEG_DTU$Age_Group == '[4-12-18]'|
                           F_DEG_DTU$Area == 'Hippocampus' & F_DEG_DTU$Age_Group == '[4-12]' |
                           F_DEG_DTU$Area == 'Hippocampus' & F_DEG_DTU$Age_Group == '[4-18]' |
                           F_DEG_DTU$Area == 'Hippocampus' & F_DEG_DTU$Age_Group == '[12-18]',1]),
  organism = 'mmusculus',evcodes = T,
  significant = T,
  multi_query = F,sources = c('GO:BP','GO:MF','GO:CC'),
  user_threshold = 0.01,
  correction_method = 'fdr')


DEG_DTU_FAD <- DEG_DTU_FAD$result[,c("query", "term_id", "term_name", "p_value", "intersection",
                             'precision','intersection_size')]
colnames(DEG_DTU_FAD) <-  c("query", "GO.ID", "Description", "p.Val",'Genes','precision',
                        'intersection_size')
DEG_DTU_FAD$FDR <- DEG_DTU_FAD$p.Val
DEG_DTU_FAD$"Class" <- 'DEG'
DEG_DTU_FAD$'Phenotype' <- "+1"
DEG_DTU_FAD <- DEG_DTU_FAD[DEG_DTU_FAD$precision >0.03,]





DEG_DTU_FAD %>% group_by(query) %>%
  group_walk(~
               write.table(data.frame(.x[,c("GO.ID", "Description", "p.Val", "FDR",'Phenotype',"Genes")]), 
                           file = paste0("NEW_GPROFILER/FAD/DEG-DTU/", unique(.y$query), "_gem.txt"),
                           sep = "\t", quote = F, row.names = F))




fad5x_gprofiler  <- rbind(DEG_FAD,DTU_FAD,DEG_DTU_FAD)
fad5x_gprofiler <- fad5x_gprofiler[fad5x_gprofiler$intersection_size > 3,]

write.csv(fad5x_gprofiler,'NEW_GPROFILER/FAD/fad5x_gprofiler2')


### TAU PER AGE

T_DEG_DTU <- readRDS('NEW_GPROFILER/refs/T_DEG_DTU.rds')

TAU_DEG <- T_DEG_DTU[T_DEG_DTU$Class == 'DEG',]

DEG_TAU <- gost(query = list(
  'HIP[4]' = TAU_DEG[TAU_DEG$Age_Group == '[4]',1],
  'HIP[17]'= TAU_DEG[TAU_DEG$Age_Group == '[17]',1]),
  organism = 'mmusculus',evcodes = T,
  significant = T,
  multi_query = F,sources = c('GO:BP','GO:MF','GO:CC'),
  user_threshold = 0.01,
  correction_method = 'fdr')



DEG_TAU <- DEG_TAU$result[,c("query", "term_id", "term_name", "p_value", "intersection",
                             'precision','intersection_size')]
colnames(DEG_TAU) <-  c("query", "GO.ID", "Description", "p.Val",'Genes','precision',
                        'intersection_size')
DEG_TAU$FDR <- DEG_TAU$p.Val
DEG_TAU$"Class" <- 'DEG'
DEG_TAU$'Phenotype' <- "+1"
DEG_TAU <- DEG_TAU[DEG_TAU$precision >0.03,]



DEG_TAU %>% group_by(query) %>%
  group_walk(~
               write.table(data.frame(.x[,c("GO.ID", "Description", "p.Val", "FDR",'Phenotype',"Genes")]), 
                           file = paste0("NEW_GPROFILER/TAU/DEG/", unique(.y$query), "_gem.txt"),
                           sep = "\t", quote = F, row.names = F))

### DTU

TAU_DTU <- T_DEG_DTU[T_DEG_DTU$Class == 'DTU',]

DTU_TAU <- gost(query = list(
  'HIP[4]' = TAU_DTU[TAU_DTU$Age_Group == '[4]',1],
  'HIP[17]'= TAU_DTU[TAU_DTU$Age_Group == '[17]',1]),
  organism = 'mmusculus',evcodes = T,
  significant = T,
  multi_query = F,sources = c('GO:BP','GO:MF','GO:CC'),
  user_threshold = 0.01,
  correction_method = 'fdr')



DTU_TAU <- DTU_TAU$result[,c("query", "term_id", "term_name", "p_value", "intersection",
                             'precision','intersection_size')]
colnames(DTU_TAU) <-  c("query", "GO.ID", "Description", "p.Val",'Genes','precision',
                        'intersection_size')
DTU_TAU$FDR <- DTU_TAU$p.Val
DTU_TAU$"Class" <- 'DTU'
DTU_TAU$'Phenotype' <- "+1"
DTU_TAU <- DTU_TAU[DTU_TAU$precision >0.03,]



DTU_TAU %>% group_by(query) %>%
  group_walk(~
               write.table(data.frame(.x[,c("GO.ID", "Description", "p.Val", "FDR",'Phenotype',"Genes")]), 
                           file = paste0("NEW_GPROFILER/TAU/DTU/", unique(.y$query), "_gem.txt"),
                           sep = "\t", quote = F, row.names = F))

### DEG-DTU


DEG_DTU_TAU <- gost(query = list(
  'HIP[4]' = T_DEG_DTU[T_DEG_DTU$Age_Group == '[4]',1],
  'HIP[17]'= T_DEG_DTU[T_DEG_DTU$Age_Group == '[17]',1]),
  organism = 'mmusculus',evcodes = T,
  significant = T,
  multi_query = F,sources = c('GO:BP','GO:MF','GO:CC'),
  user_threshold = 0.01,
  correction_method = 'fdr')



DEG_DTU_TAU <- DEG_DTU_TAU$result[,c("query", "term_id", "term_name", "p_value", "intersection",
                                     'precision','intersection_size')]
colnames(DEG_DTU_TAU) <-  c("query", "GO.ID", "Description", "p.Val",'Genes','precision',
                            'intersection_size')
DEG_DTU_TAU$FDR <- DEG_DTU_TAU$p.Val
DEG_DTU_TAU$"Class" <- 'DEG-DTU'
DEG_DTU_TAU$'Phenotype' <- "+1"
DEG_DTU_TAU <- DEG_DTU_TAU[DEG_DTU_TAU$precision >0.03,]



DEG_DTU_TAU%>% group_by(query) %>%
  group_walk(~
               write.table(data.frame(.x[,c("GO.ID", "Description", "p.Val", "FDR",'Phenotype',"Genes")]), 
                           file = paste0("NEW_GPROFILER/TAU/DEG-DTU/", unique(.y$query), "_gem.txt"),
                           sep = "\t", quote = F, row.names = F))

tau_gprofiler <- rbind(DEG_TAU,DTU_TAU,DEG_DTU_TAU)
tau_gprofiler <- tau_gprofiler[tau_gprofiler$intersection_size > 3,]

write.csv(tau_gprofiler,'NEW_GPROFILER/TAU/taud35_gprofiler.csv')

