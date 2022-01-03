library('stageR')
library('DRIMSeq')


##### Import of Isoform Switch archives

Tau_4 <- readRDS('stage_archives/refs/Mouse/age_4.fullAnalysis.rds')
Tau_4 <- Tau_4$isoformFeatures

Tau_17 <- readRDS('stage_archives/refs/Mouse/age_17.fullAnalysis.rds')
Tau_17 <- Tau_17$isoformFeatures

### Import of gene correspondece table, transcript counts and metadata

tx2gene <- read.table('stage_archives/refs/Mouse/t2g.txt',header = F)

counts_TAU <- read.csv('stage_archives/refs/Mouse/counts_TAU.csv',header = T, stringsAsFactors = F,
                       row.names = 1,check.names = F)

TAU_metadado <- read.csv('stage_archives/refs/Mouse/Tau_metadado.csv',header = T, stringsAsFactors = F,
                         row.names = 1)

### Run the stageR algorithm to mice with 4 months

TAU_metadado_4 <- TAU_metadado[TAU_metadado$age == '4_months',]

counts_TAU_4 <- counts_TAU[,TAU_metadado_4$BamFileName]
counts_TAU_4 <- counts_TAU_4[rownames(counts_TAU_4) %in% Tau_4$isoform_id,]
counts_TAU_4 <- counts_TAU_4[order(rownames(counts_TAU_4)),]

tx2gene_4 <- tx2gene[tx2gene$V1 %in% rownames(counts_TAU_4),]
tx2gene_4 <- tx2gene_4[order(tx2gene_4$V1),]

counts_TAU_4$'feature_id' <- tx2gene_4$V1
counts_TAU_4$'gene_id' <- tx2gene_4$V2

samples_TAU_4 <- data.frame(sample_id = TAU_metadado_4$BamFileName,
                            group = TAU_metadado_4$mutation)

d <- dmDSdata(counts = counts_TAU_4,samples = samples_TAU_4)

d <- dmFilter(d)
design_full <- model.matrix(~group,data=samples_TAU_4)

d <- dmPrecision(d,design = design_full)

d <- dmFit(d, design = design_full, verbose = 1)

d <- dmTest(d, coef = 'groupControl',verbose = 1)


pScreen <- d@results_gene$pvalue
names(pScreen) <- d@results_gene$gene_id

pConfirmation <-matrix(d@results_feature$pvalue)
rownames(pConfirmation) <- d@results_feature$feature_id

tx2gene_4 <- tx2gene_4[,c(1:2)]
colnames(tx2gene_4) <- c('feature_id','gene_id')

stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
                      pScreenAdjusted = FALSE, tx2gene = tx2gene_4)

stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu",
                                 alpha = 0.01)
getSignificantGenes(stageRObj)
getSignificantTx(stageRObj)

padj_4 <- getAdjustedPValues(stageRObj, order = TRUE,
                             onlySignificantGenes = T)

write.csv(padj_4,'stage_archives/results/Mouse/stage_TAU_4.csv')



### Run the stageR algorithm to mice with 17 months

TAU_metadado_17 <- TAU_metadado[TAU_metadado$age == '17_months',]

counts_TAU_17 <- counts_TAU[,TAU_metadado_17$BamFileName]
counts_TAU_17 <- counts_TAU_17[rownames(counts_TAU_17) %in% Tau_17$isoform_id,]
counts_TAU_17 <- counts_TAU_17[order(rownames(counts_TAU_17)),]

tx2gene_17 <- tx2gene[tx2gene$V1 %in% rownames(counts_TAU_17),]
tx2gene_17 <- tx2gene_17[order(tx2gene_17$V1),]

counts_TAU_17 <- counts_TAU_17[rownames(counts_TAU_17) %in% tx2gene_17$V1,]

counts_TAU_17$'feature_id' <- tx2gene_17$V1
counts_TAU_17$'gene_id' <- tx2gene_17$V2

samples_TAU_17 <- data.frame(sample_id = TAU_metadado_17$BamFileName,
                            group = TAU_metadado_17$mutation)

d <- dmDSdata(counts = counts_TAU_17,samples = samples_TAU_17)

d <- dmFilter(d)
design_full <- model.matrix(~group,data=samples_TAU_17)

d <- dmPrecision(d,design = design_full)

d <- dmFit(d, design = design_full, verbose = 1)

d <- dmTest(d, coef = 'groupControl',verbose = 1)


pScreen <- d@results_gene$pvalue
names(pScreen) <- d@results_gene$gene_id

pConfirmation <-matrix(d@results_feature$pvalue)
rownames(pConfirmation) <- d@results_feature$feature_id

tx2gene_17 <- tx2gene_17[,c(1:2)]
colnames(tx2gene_17) <- c('feature_id','gene_id')

stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
                      pScreenAdjusted = FALSE, tx2gene = tx2gene_17)

stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu",
                                 alpha = 0.01,allowNA=TRUE)
getSignificantGenes(stageRObj)
getSignificantTx(stageRObj)

padj_17 <- stageRObj@adjustedP

write.csv(padj_17,'stage_archives/results/Mouse/stage_TAU_17.csv')
