library('IsoformSwitchAnalyzeR')
library('dplyr')
library('tibble')
library('stringr')

Month_Group <- c('4','17')


## Group of Ages


ISAR_tabs <- lapply(Month_Group,function(M){

        BM <- readRDS(paste0("results/first_analysis/age_",M,"_months_aSwitchListAnalyzed.rds"))

        print('Were working')

     
	BM <- isoformSwitchAnalysisPart2(BM,
                                 pathToCPC2resultFile = paste0("results/refs_part2/cpc2/age_",M,".txt"),
                                 pathToPFAMresultFile = paste0("results/refs_part2/pfam/age_",M,"_pfam.txt"),
                                 pathToNetSurfP2resultFile = paste0("results/refs_part2/netsurfp2/age_",M,".csv"),
                                 pathToSignalPresultFile = paste0("results/refs_part2/signalp/age_",M,"_summary.signalp5"),
                                 dIFcutoff = 0.05,
                                 removeNoncodinORFs = T,
                                 outputPlots = F)

        saveRDS(BM, file = paste0("results/second_analysis/age_",M,".fullAnalysis.rds"))

        print("Its done!!!")

})

