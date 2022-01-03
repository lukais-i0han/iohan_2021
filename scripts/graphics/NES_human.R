#### Cemitoll NES graphic Human | Mouse
library('stringr')
library('ggplot2')
library('dplyr')


custom_pal <- c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                "#D6604D", "#B2182B", "#67001F")
custom_pal <- colorRampPalette(custom_pal)(200)



Nes_human_mouse <- read.csv('NEW_RESULTS/results/cemitool_01/tables/NES_human_mouse.csv',header = T,
                            stringsAsFactors = F)

Nes_human_mouse$group <- factor(Nes_human_mouse$group,
                                levels = c('AD','PSP','FAD5X','TauD35','MSBB_10','MSBB_36','CON'))
Nes_human_mouse <- Nes_human_mouse[Nes_human_mouse$padj <0.01,]

g<-ggplot(Nes_human_mouse, aes(x=modules, y=group, size=abs(NES), fill=NES)) +
  geom_point(color = "white", shape=21) +
  scale_fill_gradientn(colours=custom_pal, space = "Lab") +
  scale_size(range=c(0,30), limits=c(0, NA)) +
  guides(size="none") + facet_wrap(~condition,nrow = 2)+
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 15)) +
  scale_x_discrete(position = "top")

ALZ_TCX<- g %+% subset(Nes_human_mouse,group %in% 'AD')
PSP_TCX<- g %+% subset(Nes_human_mouse,group %in% 'PSP')
MSBB_10 <- g %+% subset(Nes_human_mouse,group %in% 'MSBB_10')
MSBB_36 <- g %+% subset(Nes_human_mouse,group %in% 'MSBB_36')
FAD5X <- g %+% subset(Nes_human_mouse,group %in% 'FAD5X')
TauD35 <- g %+% subset(Nes_human_mouse,group %in% 'TauD35')

for (i in 1:length(plot_list)) {
  ggsave(file=paste0('NEW_RESULTS/results/plot',i,'.jpeg'), plot = plot_list[[i]],
         units = 'cm' ,width = 15 ,height = 35,
         dpi =600)
  
  
}

plot_list_CER <- list(ALZ_CER,PSP_CER,PA_CER)

for (i in 1:length(plot_list_CER)) {
  ggsave(file=paste0('NEW_RESULTS/results/plot',i,'.jpeg'), plot = plot_list_CER[[i]],
         units = 'cm' ,width = 15 ,height = 35,
         dpi =1200)
  
  
}

