#### NES value to NA modules
library('CEMiTool')


ALZ_TCX <- ALZ_TCX_nes_values

cem_object <- ALZ_TCX

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

saveRDS(cem_object,'CEMtool/new_results/TCX/ALZ.rds')


### PSP


PSP <- readRDS('NEW_RESULTS/Cemtool_results/human/pval_01/PSP.rds')

cem_object <- PSP

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

saveRDS(cem_object,'CEMtool/new_results/TCX/ALZ.rds')

####  MSBB 10



MSBB_10<- readRDS('NEW_RESULTS/Cemtool_results/MSBB/ALZ_10.rds')

cem_object <- MSBB_10

z_expr <- expr_data(cem_object, filter=cem_object@parameters$filter,
                    apply_vst=T,
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
                                         scoreType = 'pos',
                                         minSize=15,
                                         maxSize=1000,
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

saveRDS(cem_object,'CEMtool/new_results/TCX/ALZ.rds')

### MSBB 36


cem_object <- readRDS('NEW_RESULTS/Cemtool_results/MSBB/ALZ_36.rds')

cem_object <- MSBB

z_expr <- expr_data(cem_object, filter=cem_object@parameters$filter,
                    apply_vst=T,
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
                                         maxSize=1000,
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

saveRDS(cem_object,'CEMtool/new_results/FAD/HIP.rds')

#### Mouse CCX

cem_object <- readRDS('CEMtool/new_results/FAD/CCX.rds')

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
                                         nPermSimple = 20000,
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

saveRDS(cem_object,'CEMtool/new_results/FAD/CCX.rds')

### Mouse HIP

cem_object <- readRDS('NEW_RESULTS/Cemtool_results/mouse_FAD5X/HIP.rds')

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
                                         maxSize=1000,
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

saveRDS(cem_object,'CEMtool/new_results/FAD/HIP.rds')


#### Taud35

Tau <- readRDS('NEW_RESULTS/Cemtool_results/mouse_TAU/Taud35.rds')
write.csv(Tau@enrichment$nes,'NEW_RESULTS/results/tables/nes_tau.csv')
write.csv(Tau@enrichment$padj,'NEW_RESULTS/results/tables/padj_tau.csv')
