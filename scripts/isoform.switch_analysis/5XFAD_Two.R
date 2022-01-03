## Isoform Part 2
library('IsoformSwitchAnalyzeR')
library('dplyr')
library('tibble')
library('stringr')

Month_Group <- c('4','12','18')


## Group of Ages


ISAR_tabs <- lapply(Month_Group,function(M){

	BM <- readRDS(paste0("results/first_analysis/Month_",M,"_aSwitchListAnalyzed.rds"))

	print('Were working')

	BM <- isoformSwitchAnalysisPart2(BM,
                                 pathToCPC2resultFile = paste0("results/refs_part2/cpc2/Month_",M,".txt"),
                                 pathToPFAMresultFile = paste0("results/refs_part2/pfam/Month_",M,"_pfam.txt"),
                                 pathToNetSurfP2resultFile = paste0("results/refs_part2/netsurfp2/Month_",M,".csv"),
                                 pathToSignalPresultFile = paste0("results/refs_part2/signal5p/Month_",M,"_summary.signalp5"),
                                 dIFcutoff = 0.05,
                                 removeNoncodinORFs = T,
                                 outputPlots = F)

	saveRDS(BM, file = paste0("results/second_analysis/Month_",M,".fullAnalysis.rds"))

	print("Its done!!!")

})

