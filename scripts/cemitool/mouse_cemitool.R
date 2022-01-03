library('CEMiTool')
library('stringr')
library('ggplot2')
library('dplyr')



### import of vsd normalised objects by DESeq2 

vsd_FAD <- readRDS('CEMtool/vsd_files/vsd_FAD.rds')
vsd_TAU <- readRDS('CEMtool/vsd_files/vsd_TAU.rds')

### import of metadado


FAD_metadado <- read.csv('CEMtool/metadado/FAD_metadado.csv',header = T, 
                         stringsAsFactors = F,row.names = 1)

Tau_metadado <- read.csv('CEMtool/metadado/Tau_metadado.csv',row.names = 1, header = T,
                         stringsAsFactors = F)

#### archives to mouse

mouse_GO <- read_gmt('CEMtool/gmt_files/Mouse_GO_AllPathways_no_GO_iea_February_05_2021_symbol.gmt')
dump_human <- str_detect(mouse_GO$term,'HUMAN')
mouse_GO <- mouse_GO[!dump_human,]

mouse_PPI <- readRDS('CEMtool/PPI/mouse_PPI.rds')


### 

vsd_FAD <- assay(vsd_FAD)

vsd_FAD <- as.data.frame(vsd_FAD)

FAD_meta <- data.frame(SampleName = FAD_metadado$SpecimenID,
                       Class = FAD_metadado$Group,
                       Region = FAD_metadado$Region)


cemitool_function <- function(vsd_object,metadado,Area){
  
  metadado <- metadado[metadado$Region == Area,]
  
  metadado <- metadado[,c(1,2)]
  
  vsd_object <- vsd_object[,colnames(vsd_object) %in% metadado$SampleName]
  
  
  cem_object <- cemitool(expr=vsd_object,
                         force_beta = T,
                         annot = metadado,
                         gsea_max_size = 2000,
                         filter_pval = 0.1,
                         merge_similar = T,
                         apply_vst = F,
                         ora_pval = 0.01,
                         filter = T)
  
  
  cem_object <- plot_gsea(cem_object)
  
  cem_object <- plot_profile(cem_object)
  
  cem_object <- mod_ora(cem_object,mouse_GO)
  cem_object <- plot_ora(cem_object)
  
  interactions_data(cem_object) <- mouse_PPI
  cem_object <- plot_interactions(cem_object)
  
  saveRDS(cem_object,file=paste0('CEMtool/new_results/FAD/',Area,'.rds'))
  write_files(cem_object,directory=paste0('CEMtool/new_results/FAD/',Area,'/'),force=T)
  
}

group_Condition <- c('CCX','HIP') ### Run for the two areas ind the data


for (i in 2) {
  
  cemitool_function(vsd_FAD,FAD_meta,group_Condition[i])
}


### TAU 


TAU_meta <- Tau_metadado[,c(1,3)]
colnames(TAU_meta) <- c('SampleName','Class')

vsd_TAU <- assay(vsd_TAU)
vsd_TAU <- as.data.frame(vsd_TAU)


cem_object <- cemitool(expr=vsd_TAU,
                         force_beta = T,
                         annot =TAU_meta,
                         gsea_max_size = 2000,
                         filter_pval = 0.1,
                         merge_similar = T,
                         apply_vst = F,
                         ora_pval = 0.01,
                         filter = T)
  
  
cem_object <- plot_gsea(cem_object)
  
cem_object <- plot_profile(cem_object)
  
cem_object <- mod_ora(cem_object,mouse_GO)
cem_object <- plot_ora(cem_object)
  
interactions_data(cem_object) <- mouse_PPI
cem_object <- plot_interactions(cem_object)
  
saveRDS(cem_object,'CEMtool/new_results/TAU/Taud35.rds')
write_files(cem_object,directory='CEMtool/new_results/TAU',force=T)
  


