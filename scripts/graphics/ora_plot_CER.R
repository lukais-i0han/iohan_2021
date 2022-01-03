##### ORA analysis CER
library('readr')
library('stringr')
library('dplyr')
library('ggplot2')

ALZ <- read_tsv('NEW_RESULTS/Cemtool_results/human/CER/ALZ/ora.tsv')

split_ontologie <- data.frame(str_split_fixed(ALZ$ID,'%',3))

ALZ[,c(2:3)] <- split_ontologie[,c(2,1)]

ALZ$Description <- str_to_lower(ALZ$Description,locale = 'en')
ALZ <- ALZ[ALZ$p.adjust < 0.01 & ALZ$Count >= 3,]

ALZ_GO <- ALZ[str_detect(ALZ$ID,'GO'),]

ALZ_GO_top10 <- ALZ_GO %>% 
  group_by(Module) %>% top_n(n=3,wt=Count)

ALZ_GO_top10 <- ALZ_GO_top10[ALZ_GO_top10$Module != 'Not.Correlated',]

ALZ_GO_top10$Module <- paste(ALZ_GO_top10$Module,'ALZ',sep = '_')

myplot_ALZ <-ggplot(ALZ_GO_top10,aes(x=Count, y=factor(Description))) +
  geom_bar(stat = 'identity')+ scale_x_continuous(breaks = c(0,15,15))+
  labs(y='Description', x=NULL) + facet_grid(~Module,shrink = T, as.table = T) +
  theme_minimal() + theme(
    axis.text.x = element_text(face = "bold", color = "black",size = 10),
    axis.title = element_text(face='bold',color='black',size=15),
    axis.text.y = element_text(face = 'bold',color = 'black',size = 10))+
  theme(legend.position = 'right',legend.justification = 'right',
        legend.text = element_text(color = "black", size = 12))


png('NEW_RESULTS/results/ALZ_CER_module_ontologie.png',height = 20,width = 25,units = 'cm',
    res = 1200)
myplot_ALZ
dev.off()
### PSP

PSP <- read_tsv('NEW_RESULTS/Cemtool_results/human/CER/PSP/ora.tsv')

split_ontologie_PSP <- data.frame(str_split_fixed(PSP$ID,'%',3))

PSP[,c(2:3)] <- split_ontologie_PSP[,c(2,1)]

PSP$Description <- str_to_lower(PSP$Description,locale = 'en')
PSP <- PSP[PSP$p.adjust < 0.01 & PSP$Count >= 3,]

PSP_GO <- PSP[str_detect(PSP$ID,'GO'),]

PSP_GO_top10 <- PSP_GO %>% 
  group_by(Module) %>% top_n(n=2,wt=Count)

PSP_GO_top10 <- PSP_GO_top10[PSP_GO_top10$Module != 'Not.Correlated',]

PSP_GO_top10$Module <- paste(PSP_GO_top10$Module,'PSP',sep = '_')

myplot_PSP <-ggplot(PSP_GO_top10,aes(x=Count, y=factor(Description))) +
  geom_bar(stat = 'identity')+ scale_x_continuous(breaks = c(0,20,10)) +
  labs(y='Description', x=NULL) + facet_grid(~Module,shrink = T, as.table = T) +
  theme_minimal() + theme(
    axis.text.x = element_text(face = "bold", color = "black",size = 10),
    axis.title = element_text(face='bold',color='black',size=15),
    axis.text.y = element_text(face = 'bold',color = 'black',size = 10))+
  theme(legend.position = 'right',legend.justification = 'right',
        legend.text = element_text(color = "black", size = 12))

png('NEW_RESULTS/results/PSP_CER_module_ontologie.png',height = 20,width = 25,units = 'cm',
    res = 1200)
myplot_PSP

dev.off()


### PA

PA <- read_tsv('NEW_RESULTS/Cemtool_results/human/CER/PA/ora.tsv')

split_ontologie_PA <- data.frame(str_split_fixed(PA$ID,'%',3))

PA[,c(2:3)] <- split_ontologie_PA[,c(2,1)]

PA$Description <- str_to_lower(PA$Description,locale = 'en')
PA <- PA[PA$p.adjust < 0.01 & PA$Count >= 3,]

PA_GO <- PA[str_detect(PA$ID,'GO'),]
PA_GO <- PA_GO[PA_GO$Module != 'Not.Correlated',]

PA_GO_top10 <- PA_GO %>% 
  group_by(Module) %>% top_n(n=3,wt=Count)


myplot_PA <-ggplot(PA_GO_top10,aes(x=Count, y=factor(Description))) +
  geom_bar(stat = 'identity')+ scale_x_continuous(breaks = c(0,10))+
  labs(y='Description', x=NULL) + facet_grid(~Module,shrink = T, as.table = T) +
  theme_minimal() + theme(
    axis.text.x = element_text(face = "bold", color = "black",size = 10),
    axis.title = element_text(face='bold',color='black',size=15),
    axis.text.y = element_text(face = 'bold',color = 'black',size = 10))+
  theme(legend.position = 'right',legend.justification = 'right',
        legend.text = element_text(color = "black", size = 12))

png('NEW_RESULTS/results/PA_CER_module_ontologie.png',height = 20,width = 20,units = 'cm',
    res = 1200)
myplot_PA

dev.off()


