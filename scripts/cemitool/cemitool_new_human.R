#### Cemitool analysis
library('CEMiTool')
library('stringr')
library('ggplot2')
library('dplyr')



### normalized matrix (vsd by DESEQ2)

vsd_mayo_TCX <- readRDS('all_data/Mayo_CER_TCX/results/new_results/TCX/vsd_mayo_TCX.rds')
vsd_mayo_CER <- readRDS('all_data/Mayo_CER_TCX/results/new_results/CER/vsd_mayo_CER.rds')



#### metadado

mayo_metadado_TCX <- read.csv('all_data/Mayo_CER_TCX/results/new_results/TCX/mayo_metadado_TCX.csv',header = T,stringsAsFactors = F,
                          row.names = 1)

mayo_metadado_CER <- read.csv('all_data/Mayo_CER_TCX/results/new_results/CER/mayo_metadado_CER.csv',header = T,stringsAsFactors = F,
              row.names = 1)

### Human

### TCX
vsd_mayo_TCX <- assay(vsd_mayo_TCX)
vsd_mayo_TCX <- as.data.frame(vsd_mayo_TCX)

mayo_meta_TCX <- data.frame(SampleName = mayo_metadado_TCX$ID,
                        Class = mayo_metadado_TCX$condition)

mayo_meta_TCX <- mayo_meta_TCX %>% mutate(Class = case_when(
  Class == 'AD' ~ 'ALZ',
  Class == 'Control' ~ 'CON',
  Class == 'Pathologic Aging' ~ 'PA',
  TRUE ~ 'PSP'
))

mayo_metadado_TCX$Class <- factor(mayo_meta_TCX$Class,levels = c('ALZ','PSP','PA','CON'))

### CER

vsd_mayo_CER <- assay(vsd_mayo_CER)
vsd_mayo_CER <- as.data.frame(vsd_mayo_CER)

mayo_meta_CER <- data.frame(SampleName = mayo_metadado_CER$ID,
                        Class = mayo_metadado_CER$condition)

mayo_meta_CER <- mayo_meta_CER %>% mutate(Class = case_when(
  Class == 'AD' ~ 'ALZ',
  Class == 'Control' ~ 'CON',
  Class == 'Pathologic Aging' ~ 'PA',
  TRUE ~ 'PSP'
))

mayo_meta_CER$Class <- factor(mayo_meta_CER$Class,levels = c('ALZ','PSP','PA','CON'))




### Human Pathway and network

human_GO<-read_gmt('CEMtool/gmt_files/Human_GO_AllPathways_no_GO_iea_February_05_2021_symbol.gmt')

int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)


### Create cem archive TCX


cemitool_function <- function(group_Condtition,vsd_object,metadado) {
  
  metadado <- metadado[metadado$Class == group_Condtition | 
                       metadado$Class == 'CON',]
  
  
  vsd_object <- vsd_object[,colnames(vsd_object) %in% metadado$SampleName]
  
  cem_object <- cemitool(expr=vsd_object,force_beta = T,
                         annot = metadado,
                         gsea_max_size = 2000,
                         filter_pval = 0.1,
                         interactions = int_df,
                         merge_similar = T,
                         apply_vst = F,
                         ora_pval = 0.01,
                         filter = T)
  
  
  cem_object <- plot_gsea(cem_object)
  
  cem_object <- plot_profile(cem_object)
  
  cem_object <- mod_ora(cem_object,human_GO)
  cem_object <- plot_ora(cem_object)
  
  interactions_data(cem_object) <- int_df
  cem_object <- plot_interactions(cem_object)
  
  saveRDS(cem_object,file=paste0('CEMtool/new_results/TCX/',group_Condtition,'.rds'))
  write_files(cem_object,directory=paste0('CEMtool/new_results/TCX/',group_Condtition,'/'),force=T)
  
}

group_Condtition <- c('PSP')

for (i in 1) {
  cemitool_function(group_Condtition[i],vsd_mayo_TCX,mayo_meta_TCX)
  
}

### Create cemtool archive

cemitool_function <- function(group_Condtition,vsd_object,metadado) {
  
  metadado <- metadado[metadado$Class == group_Condtition | 
                         metadado$Class == 'CON',]
  
  
  vsd_object <- vsd_object[,colnames(vsd_object) %in% metadado$SampleName]
  
  cem_object <- cemitool(expr=vsd_object,force_beta = T,
                         annot = metadado,
                         gsea_max_size = 2000,
                         filter_pval = 0.05,
                         interactions = int_df,
                         merge_similar = T,
                         apply_vst = F,
                         ora_pval = 0.01,
                         filter = T)
  
  
  cem_object <- plot_gsea(cem_object)
  
  cem_object <- plot_profile(cem_object)
  
  cem_object <- mod_ora(cem_object,human_GO)
  cem_object <- plot_ora(cem_object)
  
  interactions_data(cem_object) <- int_df
  cem_object <- plot_interactions(cem_object)
  
  saveRDS(cem_object,file=paste0('CEMtool/new_results/CER/',group_Condtition,'.rds'))
  write_files(cem_object,directory=paste0('CEMtool/new_results/CER/',group_Condtition,'/'),force=T)
  
}

group_Condtition <- c('ALZ','PA','PSP')

for (i in 1:3) {
  cemitool_function(group_Condtition[i],vsd_mayo_CER,mayo_meta_CER)
  
}

