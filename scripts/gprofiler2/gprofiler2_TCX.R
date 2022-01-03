### gprofiler tables too Ctyoscape
library('gprofiler2')
library('dplyr')
library('stringr')

mayo_DEG_DTU_TCX <- read.csv('NEW_GPROFILER/refs/H_DEG_DTU_TCX.csv',header = T,stringsAsFactors = F,
                         row.names = 1) ### table with genes classified by disease, expression (DEG and gDTU) and groups (A, B and C)


mayo_DEG_TCX <- mayo_DEG_DTU_TCX[mayo_DEG_DTU_TCX$Class == 'DEG',]

mayo_DTU_TCX <- mayo_DEG_DTU_TCX[mayo_DEG_DTU_TCX$Class == 'DTU',]

mayo_DEG <- gost(query = list(
  'ALZ_A' = mayo_DEG_TCX[mayo_DEG_TCX$Group == 'ALZ' & mayo_DEG_TCX$Age_Group == '[A]' | 
                         mayo_DEG_TCX$Group == 'ALZ' & mayo_DEG_TCX$Age_Group == '[AB]'|
                         mayo_DEG_TCX$Group == 'ALZ' & mayo_DEG_TCX$Age_Group == '[AC]',1],
  'ALZ_B' = mayo_DEG_TCX[mayo_DEG_TCX$Group == 'ALZ' & mayo_DEG_TCX$Age_Group == '[B]' | 
                           mayo_DEG_TCX$Group == 'ALZ' & mayo_DEG_TCX$Age_Group == '[AB]'|
                           mayo_DEG_TCX$Group == 'ALZ' & mayo_DEG_TCX$Age_Group == '[BC]',1],
  'ALZ_C' = mayo_DEG_TCX[mayo_DEG_TCX$Group == 'ALZ' & mayo_DEG_TCX$Age_Group == '[C]' | 
                           mayo_DEG_TCX$Group == 'ALZ' & mayo_DEG_TCX$Age_Group == '[AC]'|
                           mayo_DEG_TCX$Group == 'ALZ' & mayo_DEG_TCX$Age_Group == '[BC]',1],
  'ALZ[AB.AC.BC]' = mayo_DEG_TCX[mayo_DEG_TCX$Group == 'ALZ' & mayo_DEG_TCX$Age_Group == '[AB]' | 
                             mayo_DEG_TCX$Group == 'ALZ' & mayo_DEG_TCX$Age_Group == '[AC]'|
                             mayo_DEG_TCX$Group == 'ALZ' & mayo_DEG_TCX$Age_Group == '[BC]' |
                             mayo_DEG_TCX$Group == 'ALZ' & mayo_DEG_TCX$Age_Group == '[ABC]',1],
  'PSP_A' = mayo_DEG_TCX[mayo_DEG_TCX$Group == 'PSP' & mayo_DEG_TCX$Age_Group == '[A]',1],
  'PSP_B' = mayo_DEG_TCX[mayo_DEG_TCX$Group == 'PSP' & mayo_DEG_TCX$Age_Group == '[B]',1],
  'PSP_AB' = mayo_DEG_TCX[mayo_DEG_TCX$Group == 'PSP' & mayo_DEG_TCX$Age_Group == '[AB]',1],
  'PA_A' = mayo_DEG_TCX[mayo_DEG_TCX$Group == 'PA' & mayo_DEG_TCX$Age_Group == '[A]' | 
                           mayo_DEG_TCX$Group == 'PA' & mayo_DEG_TCX$Age_Group == '[AB]'|
                           mayo_DEG_TCX$Group == 'PA' & mayo_DEG_TCX$Age_Group == '[AC]',1],
  'PA_B' = mayo_DEG_TCX[mayo_DEG_TCX$Group == 'PA' & mayo_DEG_TCX$Age_Group == '[B]' | 
                           mayo_DEG_TCX$Group == 'PA' & mayo_DEG_TCX$Age_Group == '[AB]'|
                           mayo_DEG_TCX$Group == 'PA' & mayo_DEG_TCX$Age_Group == '[BC]',1],
  'PA_C' = mayo_DEG_TCX[mayo_DEG_TCX$Group == 'PA' & mayo_DEG_TCX$Age_Group == '[C]' | 
                           mayo_DEG_TCX$Group == 'PA' & mayo_DEG_TCX$Age_Group == '[AC]'|
                           mayo_DEG_TCX$Group == 'PA' & mayo_DEG_TCX$Age_Group == '[BC]',1],
  'PA[AB.AC.BC]' = mayo_DEG_TCX[mayo_DEG_TCX$Group == 'PA' & mayo_DEG_TCX$Age_Group == '[AB]' | 
                                 mayo_DEG_TCX$Group == 'PA' & mayo_DEG_TCX$Age_Group == '[AC]'|
                                 mayo_DEG_TCX$Group == 'PA' & mayo_DEG_TCX$Age_Group == '[BC]' |
                                 mayo_DEG_TCX$Group == 'PA' & mayo_DEG_TCX$Age_Group == '[ABC]',1]),
  organism = 'hsapiens',
  evcodes = T,
  significant = T,
  multi_query = F,sources = c('GO:BP','GO:MF','GO:CC'),
  user_threshold = 0.01,
  correction_method = 'fdr')



mayo_DEG <- mayo_DEG$result[,c("query", "term_id", "term_name", "p_value", "intersection",
                             'precision','intersection_size')]
colnames(mayo_DEG) <- c("query", "GO.ID", "Description", "p.Val",'Genes','precision',
                       'intersection_size')
mayo_DEG$FDR <- mayo_DEG$p.Val
mayo_DEG$"Class" <- 'DEG'
mayo_DEG$"Phenotype" <- "+1"
mayo_DEG <- mayo_DEG[mayo_DEG$precision >0.03,]
mayo_DEG <- mayo_DEG[mayo_DEG$intersection_size > 3,]

mayo_DEG %>% group_by(query) %>%
  group_walk(~
               write.table(data.frame(.x[,c("GO.ID", "Description", "p.Val", "precision",'Phenotype',"Genes")]),
                           file = paste0("NEW_GPROFILER/TCX/DEG/", unique(.y$query), "_gem.txt"),
                           sep = "\t", quote = F, row.names = F))


#### DTU

mayo_DTU <- gost(query = list(
  'ALZ_A' = mayo_DTU_TCX[mayo_DTU_TCX$Group == 'ALZ' & mayo_DTU_TCX$Age_Group == '[A]' | 
                           mayo_DTU_TCX$Group == 'ALZ' & mayo_DTU_TCX$Age_Group == '[AB]'|
                           mayo_DTU_TCX$Group == 'ALZ' & mayo_DTU_TCX$Age_Group == '[AC]',1],
  'ALZ_B' = mayo_DTU_TCX[mayo_DTU_TCX$Group == 'ALZ' & mayo_DTU_TCX$Age_Group == '[B]' | 
                           mayo_DTU_TCX$Group == 'ALZ' & mayo_DTU_TCX$Age_Group == '[AB]'|
                           mayo_DTU_TCX$Group == 'ALZ' & mayo_DTU_TCX$Age_Group == '[BC]',1],
  'ALZ_C' = mayo_DTU_TCX[mayo_DTU_TCX$Group == 'ALZ' & mayo_DTU_TCX$Age_Group == '[C]' | 
                           mayo_DTU_TCX$Group == 'ALZ' & mayo_DTU_TCX$Age_Group == '[AC]'|
                           mayo_DTU_TCX$Group == 'ALZ' & mayo_DTU_TCX$Age_Group == '[BC]',1],
  'ALZ[AB.AC.BC]' = mayo_DTU_TCX[mayo_DTU_TCX$Group == 'ALZ' & mayo_DTU_TCX$Age_Group == '[AB]' | 
                                   mayo_DTU_TCX$Group == 'ALZ' & mayo_DTU_TCX$Age_Group == '[AC]'|
                                   mayo_DTU_TCX$Group == 'ALZ' & mayo_DTU_TCX$Age_Group == '[BC]' |
                                   mayo_DTU_TCX$Group == 'ALZ' & mayo_DTU_TCX$Age_Group == '[ABC]',1],
  'PSP_A' = mayo_DTU_TCX[mayo_DTU_TCX$Group == 'PSP' & mayo_DTU_TCX$Age_Group == '[A]',1],
  'PSP_B' = mayo_DTU_TCX[mayo_DTU_TCX$Group == 'PSP' & mayo_DTU_TCX$Age_Group == '[B]',1],
  'PSP_AB' = mayo_DTU_TCX[mayo_DTU_TCX$Group == 'PSP' & mayo_DTU_TCX$Age_Group == '[AB]',1],
  'PA_A' = mayo_DTU_TCX[mayo_DTU_TCX$Group == 'PA' & mayo_DTU_TCX$Age_Group == '[A]' | 
                          mayo_DTU_TCX$Group == 'PA' & mayo_DTU_TCX$Age_Group == '[AB]'|
                          mayo_DTU_TCX$Group == 'PA' & mayo_DTU_TCX$Age_Group == '[AC]',1],
  'PA_B' = mayo_DTU_TCX[mayo_DTU_TCX$Group == 'PA' & mayo_DTU_TCX$Age_Group == '[B]' | 
                          mayo_DTU_TCX$Group == 'PA' & mayo_DTU_TCX$Age_Group == '[AB]'|
                          mayo_DTU_TCX$Group == 'PA' & mayo_DTU_TCX$Age_Group == '[BC]',1],
  'PA_C' = mayo_DTU_TCX[mayo_DTU_TCX$Group == 'PA' & mayo_DTU_TCX$Age_Group == '[C]' | 
                          mayo_DTU_TCX$Group == 'PA' & mayo_DTU_TCX$Age_Group == '[AC]'|
                          mayo_DTU_TCX$Group == 'PA' & mayo_DTU_TCX$Age_Group == '[BC]',1],
  'PA[AB.AC.BC]' = mayo_DTU_TCX[mayo_DTU_TCX$Group == 'PA' & mayo_DTU_TCX$Age_Group == '[AB]' | 
                                  mayo_DTU_TCX$Group == 'PA' & mayo_DTU_TCX$Age_Group == '[AC]'|
                                  mayo_DTU_TCX$Group == 'PA' & mayo_DTU_TCX$Age_Group == '[BC]' |
                                  mayo_DTU_TCX$Group == 'PA' & mayo_DTU_TCX$Age_Group == '[ABC]',1]),
  organism = 'hsapiens',
  evcodes = T,
  significant = T,
  multi_query = F,sources = c('GO:BP','GO:MF','GO:CC'),
  user_threshold = 0.01,
  correction_method = 'fdr')




mayo_DTU <- mayo_DTU$result[,c("query", "term_id", "term_name", "p_value", "intersection",
                                   'precision','intersection_size')]
colnames(mayo_DTU) <- c("query", "GO.ID", "Description", "p.Val",'Genes','precision',
                            'intersection_size')
mayo_DTU$FDR <- mayo_DTU$p.Val
mayo_DTU$"Class" <- 'DTU'
mayo_DTU$"Phenotype" <- "+1"
mayo_DTU <- mayo_DTU[mayo_DTU$precision >0.03,]


mayo_DTU %>% group_by(query) %>%
  group_walk(~
               write.table(data.frame(.x[,c("GO.ID", "Description", "p.Val", "precision",'Phenotype',"Genes")]),
                           file = paste0("NEW_GPROFILER/TCX/DTU/", unique(.y$query), "_gem.txt"),
                           sep = "\t", quote = F, row.names = F))

### DEG-DTU


mayo_DEG_DTU <- gost(query = list(
  'ALZ_A' = mayo_DEG_DTU_TCX[mayo_DEG_DTU_TCX$Group == 'ALZ' & mayo_DEG_DTU_TCX$Age_Group == '[A]' | 
                           mayo_DEG_DTU_TCX$Group == 'ALZ' & mayo_DEG_DTU_TCX$Age_Group == '[AB]'|
                           mayo_DEG_DTU_TCX$Group == 'ALZ' & mayo_DEG_DTU_TCX$Age_Group == '[AC]',1],
  'ALZ_B' = mayo_DEG_DTU_TCX[mayo_DEG_DTU_TCX$Group == 'ALZ' & mayo_DEG_DTU_TCX$Age_Group == '[B]' | 
                           mayo_DEG_DTU_TCX$Group == 'ALZ' & mayo_DEG_DTU_TCX$Age_Group == '[AB]'|
                           mayo_DEG_DTU_TCX$Group == 'ALZ' & mayo_DEG_DTU_TCX$Age_Group == '[BC]',1],
  'ALZ_C' = mayo_DEG_DTU_TCX[mayo_DEG_DTU_TCX$Group == 'ALZ' & mayo_DEG_DTU_TCX$Age_Group == '[C]' | 
                           mayo_DEG_DTU_TCX$Group == 'ALZ' & mayo_DEG_DTU_TCX$Age_Group == '[AC]'|
                           mayo_DEG_DTU_TCX$Group == 'ALZ' & mayo_DEG_DTU_TCX$Age_Group == '[BC]',1],
  'ALZ[AB.AC.BC]' = mayo_DEG_DTU_TCX[mayo_DEG_DTU_TCX$Group == 'ALZ' & mayo_DEG_DTU_TCX$Age_Group == '[AB]' | 
                                   mayo_DEG_DTU_TCX$Group == 'ALZ' & mayo_DEG_DTU_TCX$Age_Group == '[AC]'|
                                   mayo_DEG_DTU_TCX$Group == 'ALZ' & mayo_DEG_DTU_TCX$Age_Group == '[BC]' |
                                   mayo_DEG_DTU_TCX$Group == 'ALZ' & mayo_DEG_DTU_TCX$Age_Group == '[ABC]',1],
  'PSP_A' = mayo_DEG_DTU_TCX[mayo_DEG_DTU_TCX$Group == 'PSP' & mayo_DEG_DTU_TCX$Age_Group == '[A]',1],
  'PSP_B' = mayo_DEG_DTU_TCX[mayo_DEG_DTU_TCX$Group == 'PSP' & mayo_DEG_DTU_TCX$Age_Group == '[B]',1],
  'PSP_AB' = mayo_DEG_DTU_TCX[mayo_DEG_DTU_TCX$Group == 'PSP' & mayo_DEG_DTU_TCX$Age_Group == '[AB]',1],
  'PA_A' = mayo_DEG_DTU_TCX[mayo_DEG_DTU_TCX$Group == 'PA' & mayo_DEG_DTU_TCX$Age_Group == '[A]' | 
                          mayo_DEG_DTU_TCX$Group == 'PA' & mayo_DEG_DTU_TCX$Age_Group == '[AB]'|
                          mayo_DEG_DTU_TCX$Group == 'PA' & mayo_DEG_DTU_TCX$Age_Group == '[AC]',1],
  'PA_B' = mayo_DEG_DTU_TCX[mayo_DEG_DTU_TCX$Group == 'PA' & mayo_DEG_DTU_TCX$Age_Group == '[B]' | 
                          mayo_DEG_DTU_TCX$Group == 'PA' & mayo_DEG_DTU_TCX$Age_Group == '[AB]'|
                          mayo_DEG_DTU_TCX$Group == 'PA' & mayo_DEG_DTU_TCX$Age_Group == '[BC]',1],
  'PA_C' = mayo_DEG_DTU_TCX[mayo_DEG_DTU_TCX$Group == 'PA' & mayo_DEG_DTU_TCX$Age_Group == '[C]' | 
                          mayo_DEG_DTU_TCX$Group == 'PA' & mayo_DEG_DTU_TCX$Age_Group == '[AC]'|
                          mayo_DEG_DTU_TCX$Group == 'PA' & mayo_DEG_DTU_TCX$Age_Group == '[BC]',1],
  'PA[AB.AC.BC]' = mayo_DEG_DTU_TCX[mayo_DEG_DTU_TCX$Group == 'PA' & mayo_DEG_DTU_TCX$Age_Group == '[AB]' | 
                                  mayo_DEG_DTU_TCX$Group == 'PA' & mayo_DEG_DTU_TCX$Age_Group == '[AC]'|
                                  mayo_DEG_DTU_TCX$Group == 'PA' & mayo_DEG_DTU_TCX$Age_Group == '[BC]' |
                                  mayo_DEG_DTU_TCX$Group == 'PA' & mayo_DEG_DTU_TCX$Age_Group == '[ABC]',1]),
  organism = 'hsapiens',
  evcodes = T,
  significant = T,
  multi_query = F,sources = c('GO:BP','GO:MF','GO:CC'),
  user_threshold = 0.01,
  correction_method = 'fdr')



mayo_DEG_DTU <- mayo_DEG_DTU$result[,c("query", "term_id", "term_name", "p_value", "intersection",
                                   'precision','intersection_size')]
colnames(mayo_DEG_DTU) <- c("query", "GO.ID", "Description", "p.Val",'Genes','precision',
                            'intersection_size')
mayo_DEG_DTU$FDR <- mayo_DEG_DTU$p.Val
mayo_DEG_DTU$"Class" <- 'DEG-DTU'
mayo_DEG_DTU$"Phenotype" <- "+1"
mayo_DEG_DTU <- mayo_DEG_DTU[mayo_DEG_DTU$precision >0.03,]


mayo_DEG_DTU %>% group_by(query) %>%
  group_walk(~
               write.table(data.frame(.x[,c("GO.ID", "Description", "p.Val", "precision",'Phenotype',"Genes")]),
                           file = paste0("NEW_GPROFILER/TCX/DEG-DTU/", unique(.y$query), "_gem.txt"),
                           sep = "\t", quote = F, row.names = F))


mayo_gprofiler <- rbind(mayo_DEG,mayo_DTU,mayo_DEG_DTU)
mayo_gprofiler <- mayo_gprofiler[mayo_gprofiler$intersection_size > 3,]
write.csv(mayo_gprofiler,'NEW_GPROFILER/TCX/mayo_gprofiler_TCX.csv')




