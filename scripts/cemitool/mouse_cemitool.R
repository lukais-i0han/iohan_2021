library('CEMiTool')
library('stringr')
library('ggplot2')
library('dplyr')



### import of vsd

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

group_Condition <- c('HIP')


for (i in 1) {
  
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
  

#### Hippocampus. Run the gsea part mannually

metadado <- FAD_metadado[FAD_metadado$Region == 'HIP',]

metadado <- metadado[,c(1,2)]

vsd_object <- vsd_FAD[,colnames(vsd_FAD) %in% FAD_meta_HIP$SampleName]

cem_object <- cemitool(expr=vsd_object,
                       force_beta = T,
                       annot = FAD_meta_HIP,
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

z_expr <- expr_data(cem_object, filter=cem_object@parameters$filter,
                    apply_vst=F,
                    filter_pval=cem_object@parameters$filter_pval)


modules <- unique(cem_object@module[, 'modules'])
gene_sets <- lapply(modules, function(mod){
  return(cem_object@module[cem_object@module[, 'modules']==mod, 'genes'])
})
names(gene_sets) <- modules


annot <- cem_object@sample_annotation
class_col <- cem_object@class_column
sample_col <- cem_object@sample_name_column
classes <- unique(annot[, class_col])


gsea_list <- lapply(classes, function(class_group){
  
  # samples of class == class_group
  class_samples <- annot[annot[, class_col]==class_group, sample_col]
  
  # genes ranked by rank_method
  genes_ranked <- apply(z_expr[, class_samples, drop=FALSE], 1, tolower('mean'))
  genes_ranked <- sort(genes_ranked, decreasing=TRUE)
  
  # BiocParallel setting up
  BiocParallel::register(BiocParallel::SerialParam())
  
  gsea_results <- fgsea::fgseaMultilevel(pathways=gene_sets,
                               stats=genes_ranked,
                               minSize=15,
                               maxSize=2000,
                               nPermSimple = 10000,
                               eps=0,
                               nproc=0)
  data.table::setDF(gsea_results)
  gsea_results[, 'leadingEdge'] <- unlist(lapply(gsea_results[, 'leadingEdge'],
                                                 function(ledges){
                                                   ledges <- paste(ledges, collapse=",")
                                                 }))
  columns <- colnames(gsea_results)
  colnames(gsea_results) <- c(columns[1], paste0(columns[-1], "_", class_group))
  return(gsea_results)
})


all_classes_df <- Reduce(function(x,y) {
  merge(x,y, all=TRUE, by='pathway')
}, gsea_list)

# separating ES / NES / pval
patterns <- list('es'='^ES_','nes'='^NES_', 'padj'='^padj_')
out_gsea <- lapply(patterns, function(pattern) {
  desired_stat <- all_classes_df[, c('pathway',
                                     grep(pattern, colnames(all_classes_df),value=TRUE))]
  colnames(desired_stat) <- gsub(pattern, '', colnames(desired_stat))
  return(desired_stat)
})

names(out_gsea) <- names(patterns)

cem_object@enrichment <- out_gsea

cem_object <- plot_gsea(cem_object)

cem_object <- plot_profile(cem_object)

cem_object <- mod_ora(cem_object,mouse_GO)
cem_object <- plot_ora(cem_object)

interactions_data(cem_object) <- mouse_PPI
cem_object <- plot_interactions(cem_object)

saveRDS(cem_object,file=paste0('CEMtool/new_results/FAD/HIP.rds'))
write_files(cem_object,directory=paste0('CEMtool/new_results/FAD/HIP'),force=T)

###

#### Cortex. Run the gsea part mannually

metadado <- FAD_metadado[FAD_metadado$Region == 'CCX',]

metadado <- metadado[,c(9,6)]
colnames(metadado) <- c('SampleName','Class')
vsd_object <- vsd_FAD[,colnames(vsd_FAD) %in% metadado$SampleName]

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

z_expr <- expr_data(cem_object, filter=cem_object@parameters$filter,
                    apply_vst=F,
                    filter_pval=cem_object@parameters$filter_pval)


modules <- unique(cem_object@module[, 'modules'])
gene_sets <- lapply(modules, function(mod){
  return(cem_object@module[cem_object@module[, 'modules']==mod, 'genes'])
})
names(gene_sets) <- modules


annot <- cem_object@sample_annotation
class_col <- cem_object@class_column
sample_col <- cem_object@sample_name_column
classes <- unique(annot[, class_col])


gsea_list <- lapply(classes, function(class_group){
  
  # samples of class == class_group
  class_samples <- annot[annot[, class_col]==class_group, sample_col]
  
  # genes ranked by rank_method
  genes_ranked <- apply(z_expr[, class_samples, drop=FALSE], 1, tolower('mean'))
  genes_ranked <- sort(genes_ranked, decreasing=TRUE)
  
  # BiocParallel setting up
  BiocParallel::register(BiocParallel::SerialParam())
  
  gsea_results <- fgsea::fgseaMultilevel(pathways=gene_sets,
                                         stats=genes_ranked,
                                         minSize=15,
                                         maxSize=2000,
                                         nPermSimple = 10000,
                                         eps=0,
                                         nproc=0)
  data.table::setDF(gsea_results)
  gsea_results[, 'leadingEdge'] <- unlist(lapply(gsea_results[, 'leadingEdge'],
                                                 function(ledges){
                                                   ledges <- paste(ledges, collapse=",")
                                                 }))
  columns <- colnames(gsea_results)
  colnames(gsea_results) <- c(columns[1], paste0(columns[-1], "_", class_group))
  return(gsea_results)
})


all_classes_df <- Reduce(function(x,y) {
  merge(x,y, all=TRUE, by='pathway')
}, gsea_list)

# separating ES / NES / pval
patterns <- list('es'='^ES_','nes'='^NES_', 'padj'='^padj_')
out_gsea <- lapply(patterns, function(pattern) {
  desired_stat <- all_classes_df[, c('pathway',
                                     grep(pattern, colnames(all_classes_df),value=TRUE))]
  colnames(desired_stat) <- gsub(pattern, '', colnames(desired_stat))
  return(desired_stat)
})

names(out_gsea) <- names(patterns)

cem_object@enrichment <- out_gsea

saveRDS(cem_object,file=paste0('CEMtool/new_results/FAD/CCX.rds'))
write_files(cem_object,directory=paste0('CEMtool/new_results/FAD/CCX'),force=T)
