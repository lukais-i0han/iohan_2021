##### ORA analys
library('readr')
library('stringr')
library('dplyr')
library('ggplot2')
library('scales')

ALZ <- readRDS('NEW_RESULTS/Cemtool_results/human/pval_01/ALZ.rds')

ALZ_ora <- ALZ@ora

ALZ_ora <- ALZ_ora[str_detect(ALZ_ora$ID,'GO:'),]

split_ontologie_ALZ <- data.frame(str_split_fixed(ALZ_ora$ID,'%',3))

ALZ_ora[,c(2:3)] <- split_ontologie_ALZ[,c(2,1)]

ALZ_ora$Description <- str_to_lower(ALZ_ora$Description,locale = 'en')

ALZ_ora <- ALZ_ora[ALZ_ora$Module != 'Not.Correlated',]
ALZ_ora <- ALZ_ora[ALZ_ora$p.adjust <0.01,]
write.csv(ALZ_ora,'NEW_RESULTS/results/cemitool_01/ora_plots/AD_ora.csv')


ALZ_ora <- ALZ_ora %>% group_by(Module) %>%
  slice_min(n = 5,order_by = p.adjust)

ALZ_ora$Description <- factor(ALZ_ora$Description,
                              levels = rev(levels(factor(ALZ_ora$Description))))

myplot_ALZ <-ggplot(ALZ_ora,aes(x=-log10(p.adjust), y=Description,fill=Module)) +
  geom_bar(stat = 'identity',alpha=1, width=0.8)+
  labs(y='Description', x=expression(-log[10](FDR)))  + 
  scale_fill_manual(values=c('#F8766D','#CD9600','#7CAE00',
                              '#00BE67','#00BFC4','#00A9FF',
                              '#C77CFF','#FF61CC'))+
  theme_classic2() + theme(
    axis.text.x = element_text(face = "bold", color = "black",size = 10),
    axis.title = element_text(face='bold',color='black',size=12),
    axis.text.y = element_text(face = 'bold',color = 'black',size = 10))+
  theme(legend.position = 'right',legend.justification = 'top',
        legend.text = element_text(color = "black", size = 12))




png('NEW_RESULTS/results/cemitool_01/mayo_ALZ_modules_top.png',height = 20,width = 20,units = 'cm',
    res = 1200)
myplot_ALZ
dev.off()
### PSP

PSP <- readRDS('NEW_RESULTS/Cemtool_results/human/pval_01/PSP.rds')

PSP_ora <- PSP@ora

PSP_ora <- PSP_ora[str_detect(PSP_ora$ID,'GO:'),]

split_ontologie_PSP <- data.frame(str_split_fixed(PSP_ora$ID,'%',3))

PSP_ora[,c(2:3)] <- split_ontologie_PSP[,c(2,1)]

PSP_ora$Description <- str_to_lower(PSP_ora$Description,locale = 'en')
PSP_ora <- PSP_ora[PSP_ora$Module != 'Not.Correlated',]
PSP_ora <- PSP_ora[PSP_ora$p.adjust < 0.01,]

PSP_ora$Description <- factor(PSP_ora$Description,
        levels = rev(levels(factor(PSP_ora$Description))))

write.csv(PSP_ora,'NEW_RESULTS/results/cemitool_01/ora_plots/PSP_ora.csv')

PSP_ora_top <- PSP_ora %>% group_by(Module) %>%
  slice_min(order_by = p.adjust,n=5)

myplot_PSP <-ggplot(PSP_ora_top,aes(x=-log10(p.adjust), y=Description,fill=Module)) +
  geom_bar(stat = 'identity',alpha=1, width=0.8)+
  scale_fill_manual(values=c('#F8766D','#CD9600','#7CAE00',
                             '#00BE67','#00BFC4','#00A9FF',
                             '#C77CFF','#FF61CC'))+
  labs(y='Description', x=expression(-log[10](FDR)))  + 
  theme_classic2() + theme(
    axis.text.x = element_text(face = "bold", color = "black",size = 10),
    axis.title = element_text(face='bold',color='black',size=12),
    axis.text.y = element_text(face = 'bold',color = 'black',size = 10))+
  theme(legend.position = 'right',legend.justification = 'top',
        legend.text = element_text(color = "black", size = 12))




png('NEW_RESULTS/results/cemitool_01/mayo_PSP_modules_top.png',height = 20,width = 20,units = 'cm',
    res = 1200)
myplot_PSP
dev.off()

#### MSBB


MSBB <- readRDS('NEW_RESULTS/Cemtool_results/MSBB/ALZ_36.rds')
write.csv(MSBB@module,'NEW_RESULTS/results/cemitool_01/modules_MSBB_36.csv')

MSBB_ora <- MSBB@ora

MSBB_ora <- MSBB_ora[str_detect(MSBB_ora$ID,'GO:'),]

split_ontologie_MSBB <- data.frame(str_split_fixed(MSBB_ora$ID,'%',3))

MSBB_ora[,c(2:3)] <- split_ontologie_MSBB[,c(2,1)]

MSBB_ora$Description <- str_to_lower(MSBB_ora$Description,locale = 'en')
MSBB_ora <- MSBB_ora[MSBB_ora$Module != 'Not.Correlated',]
MSBB_ora <- MSBB_ora[MSBB_ora$p.adjust < 0.01,]

MSBB_ora$Description <- factor(MSBB_ora$Description,
                              levels = rev(levels(factor(MSBB_ora$Description))))

write.csv(MSBB_ora, "NEW_RESULTS/results/cemitool_01/MSBB_ora_36.csv")

MSBB_ora_top <- MSBB_ora %>% group_by(Module) %>%
  slice_min(order_by = p.adjust,n=10)

myplot_MSBB <-ggplot(MSBB_ora_top,aes(x=-log10(p.adjust), y=Description,fill=Module)) +
  geom_bar(stat = 'identity',alpha=1, width=0.8)+
  labs(y='Description', x=expression(-log[10](FDR)))  + 
  scale_fill_manual(values=c('#F8766D','#CD9600','#7CAE00',
                             '#00BE67','#00BFC4','#00A9FF','#FF61CC'))+
  theme_classic2() + theme(
    axis.text.x = element_text(face = "bold", color = "black",size = 10),
    axis.title = element_text(face='bold',color='black',size=12),
    axis.text.y = element_text(face = 'bold',color = 'black',size = 10))+
  theme(legend.position = 'right',legend.justification = 'top',
        legend.text = element_text(color = "black", size = 12))




png('NEW_RESULTS/results/cemitool_01/MSBB_36_modules_top.png',
    height = 20,width = 30,units = 'cm',
    res = 1200)
myplot_MSBB
dev.off()
write.csv(MSBB_ora,'NEW_RESULTS/results/cemitool_01/ora_plots/MSBB_36_ora.csv')

### PA

PA <- read_tsv('NEW_RESULTS/Cemtool_results/human/TCX/PA/ora.tsv')

split_ontologie_PA <- data.frame(str_split_fixed(PA$ID,'%',3))

PA[,c(2:3)] <- split_ontologie_PA[,c(2,1)]

PA$Description <- str_to_lower(PA$Description,locale = 'en')
PA <- PA[PA$p.adjust < 0.01 & PA$Count >= 3,]

PA_GO <- PA[str_detect(PA$ID,'GO'),]

PA_GO_top10 <- PA_GO %>% 
  group_by(Module) %>% top_n(n=4,wt=Count)

PA_GO_top10 <- PA_GO_top10[PA_GO_top10$Module != 'Not.Correlated',]

myplot_PA <-ggplot(PA_GO_top10,aes(x=Count, y=factor(Description))) +
  geom_bar(stat = 'identity')+ scale_x_continuous(breaks = c(0,75,75))+
  labs(y='Description', x=NULL) + facet_grid(~Module,shrink = T, as.table = T) +
  theme_minimal() + theme(
    axis.text.x = element_text(face = "bold", color = "black",size = 10),
    axis.title = element_text(face='bold',color='black',size=15),
    axis.text.y = element_text(face = 'bold',color = 'black',size = 10))+
  theme(legend.position = 'right',legend.justification = 'right',
        legend.text = element_text(color = "black", size = 12))

png('NEW_RESULTS/results/PA_module_TCX_ontologie.png',height = 20,width = 20,units = 'cm',
    res = 1200)
myplot_PA

dev.off()


### FAD5X CCX

FAD5X_CCX <- read_tsv('NEW_RESULTS/Cemtool_results/mouse_FAD5X/CCX/ora.tsv')

split_ontologie_CCX<- data.frame(str_split_fixed(FAD5X_CCX$ID,'%',3))

FAD5X_CCX[,c(2:3)] <- split_ontologie_CCX[,c(2,1)]

FAD5X_CCX$Description <- str_to_lower(FAD5X_CCX$Description,locale = 'en')
FAD5X_CCX <- FAD5X_CCX[FAD5X_CCX$p.adjust < 0.01 & FAD5X_CCX$Count > 3,]

FAD5X_CCX_GO <- FAD5X_CCX[str_detect(FAD5X_CCX$ID,'GO'),]

FAD5X_GO_top10 <- FAD5X_CCX_GO %>% 
  group_by(Module) %>% top_n(n=15,wt=Count)

myplot_CCX <-ggplot(FAD5X_GO_top10,aes(x=Count, y=factor(Description))) +
  geom_bar(stat = 'identity')+ scale_x_continuous(breaks = seq(0,30,by=15))+
  labs(y='Description', x=NULL) + facet_grid(~Module,shrink = T, as.table = T) +
  theme_minimal() + theme(
    axis.text.x = element_text(face = "bold", color = "black",size = 10),
    axis.title = element_text(face='bold',color='black',size=15),
    axis.text.y = element_text(face = 'bold',color = 'black',size = 10))+
  theme(legend.position = 'right',legend.justification = 'right',
        legend.text = element_text(color = "black", size = 12))

png('NEW_RESULTS/results/CCX_module_ontologie.png',height = 20,width = 25,units = 'cm',
    res = 1200)
myplot_CCX

dev.off()

### FAD5X HIP


FAD5X_HIP <- readRDS('NEW_RESULTS/Cemtool_results/mouse_FAD5X/HIP.rds')

FAD5X_HIP_ora <- FAD5X_HIP@ora
FAD5X_HIP_ora <- FAD5X_HIP_ora[str_detect(FAD5X_HIP_ora$ID,'GO:'),]

split_ontologie_FAD5X_HIP <- data.frame(str_split_fixed(FAD5X_HIP_ora$ID,'%',3))

FAD5X_HIP_ora[,c(2:3)] <- split_ontologie_FAD5X_HIP[,c(2,1)]

FAD5X_HIP_ora$Description <- str_to_lower(FAD5X_HIP_ora$Description,locale = 'en')
FAD5X_HIP_ora <- FAD5X_HIP_ora[FAD5X_HIP_ora$Module != 'Not.Correlated',]
FAD5X_HIP_ora <- FAD5X_HIP_ora[FAD5X_HIP_ora$p.adjust < 0.01,]
write.csv(FAD5X_HIP_ora,'NEW_RESULTS/results/cemitool_01/ora_plots/FAD5X_ora.csv')


FAD5X_HIP_ora <- FAD5X_HIP_ora[FAD5X_HIP_ora$Module %in% c('M1','M5','M7','M8'),]

FAD5X_HIP_ora_top <- FAD5X_HIP_ora %>% group_by(Module) %>%
  slice_min(order_by = p.adjust,n=10)

FAD5X_HIP_ora_top$Description <- factor(FAD5X_HIP_ora_top$Description,
          levels = rev(levels(factor(FAD5X_HIP_ora_top$Description))))



myplot_FAD5X <-ggplot(FAD5X_HIP_ora_top,aes(x=-log10(p.adjust), y=Description,fill=Module)) +
  geom_bar(stat = 'identity',alpha=1, width=0.8)+
  scale_fill_manual(values=c('#F8766D',
                             '#00BFC4',
                             '#C77CFF','#FF61CC'))+
  labs(y='Description',x=expression(-log[10](FDR)))  + 
  theme_classic2() + theme(
    axis.text.x = element_text(face = "bold", color = "black",size = 10),
    axis.title = element_text(face='bold',color='black',size=12),
    axis.text.y = element_text(face = 'bold',color = 'black',size = 10))+
  theme(legend.position = 'right',legend.justification = 'top',
        legend.text = element_text(color = "black", size = 12))




png('NEW_RESULTS/results/cemitool_01/FAD5X_HIP_modules_top.png',height = 20,width = 20,units = 'cm',
    res = 1200)
myplot_FAD5X
dev.off()

### TAuD35

TAU_HIP <- readRDS('NEW_RESULTS/Cemtool_results/mouse_TAU/Taud35.rds')


TAU_HIP_ora <- TAU_HIP@ora

TAU_HIP_ora <- TAU_HIP_ora[str_detect(TAU_HIP_ora$ID,'GO:'),]

split_ontologie_TAU_HIP <- data.frame(str_split_fixed(TAU_HIP_ora$ID,'%',3))

TAU_HIP_ora[,c(2:3)] <- split_ontologie_TAU_HIP[,c(2,1)]

TAU_HIP_ora$Description <- str_to_lower(TAU_HIP_ora$Description,locale = 'en')
TAU_HIP_ora <- TAU_HIP_ora[TAU_HIP_ora$Module != 'Not.Correlated',]
TAU_HIP_ora <- TAU_HIP_ora[TAU_HIP_ora$p.adjust < 0.01,]
write.csv(TAU_HIP_ora,'NEW_RESULTS/results/cemitool_01/ora_plots/TauD35_ora.csv')


TAU_HIP_ora_top <- TAU_HIP_ora %>% group_by(Module) %>%
  slice_min(order_by = p.adjust,n=10)

TAU_HIP_ora_top$Description <- factor(TAU_HIP_ora_top$Description,
                                    levels = rev(levels(factor(TAU_HIP_ora_top$Description))))
 

myplot_TAU <-ggplot(TAU_HIP_ora_top,aes(x=-log10(p.adjust), y=Description,fill=Module)) +
  geom_bar(stat = 'identity',alpha=1, width=0.8)+
  scale_fill_manual(values=c('#F8766D','#7CAE00',
                             '#00BFC4'))+
  labs(y='Description', x=expression(-log[10](FDR)))  + 
  theme_classic2() + theme(
    axis.text.x = element_text(face = "bold", color = "black",size = 10),
    axis.title = element_text(face='bold',color='black',size=12),
    axis.text.y = element_text(face = 'bold',color = 'black',size = 10))+
  theme(legend.position = 'right',legend.justification = 'top',
        legend.text = element_text(color = "black", size = 12))




png('NEW_RESULTS/results/cemitool_01/TAU_HIP_modules_top.png',height = 20,width = 20,units = 'cm',
    res = 1200)
myplot_TAU
dev.off()


theme_classic2 <- function(base_size = 15, base_family = ""){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.border     = element_blank(),
      panel.background = element_rect(color = 'black',size = 0.5),
      axis.line        = element_line(colour = "black"),
      panel.grid.minor.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.background = element_blank(),
      legend.key       = element_blank()
      
    )
}
