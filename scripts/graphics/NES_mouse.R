#### Cemitoll NES graphic Mouse

library('CEMiTool')
library('stringr')
library('ggplot2')
library('dplyr')

FAD_CCX <- readRDS('NEW_RESULTS/Cemtool_results/mouse_FAD5X/CCX.rds')
FAD_HIP <- readRDS('NEW_RESULTS/Cemtool_results/mouse_FAD5X/HIP.rds')
TAU <- readRDS('NEW_RESULTS/Cemtool_results/mouse_TAU/Taud35.rds')

NES_FAD_CCX <- FAD_CCX@enrichment$nes
NES_FAD_CCX$'Area' <- 'CCX'
NES_FAD_CCX$'padj_CCX' <- FAD_CCX@enrichment$padj$Alzheimer
NES_FAD_CCX$'padj_CON' <- FAD_CCX@enrichment$padj$Control

NES_FAD_HIP <- FAD_HIP@enrichment$nes
colnames(NES_FAD_HIP) <- c('pathway','Alzheimer','Control')
NES_FAD_HIP$'Area' <- 'HIP'
NES_FAD_HIP$'padj_HIP' <- FAD_HIP@enrichment$padj$Alzheimer
NES_FAD_HIP$'padj_CON' <- FAD_HIP@enrichment$padj$Control

NES_TAU <- TAU@enrichment$nes
NES_TAU$'Area' <- 'HIP'
NES_TAU$'padj_TAU' <- TAU@enrichment$padj$Alzheimer
NES_TAU$'padj_CON' <- TAU@enrichment$padj$Control

NES_FAD <- list(NES_FAD_CCX,NES_FAD_HIP)

nes_tabs <- lapply(NES_FAD,function(NES) {
  
  NES <- data.frame(modules=rep(NES$pathway,2),NES=c(NES[,2],NES[,3]),
                    Group = c(rep(colnames(NES)[2],length(NES$pathway)),
                              rep(colnames(NES)[3],length(NES$pathway))),
                    Model = rep('FAD',length(NES$pathway)*2),
                    Condition =rep(colnames(NES)[2],length(NES$pathway)*2),
                    Area = NES$Area,
                    padj = c(NES[,5],NES[,6]))
})


NES_TAU <- data.frame(modules=rep(NES_TAU$pathway,2),NES=c(NES_TAU[,2],NES_TAU[,3]),
              Group = c(rep(colnames(NES_TAU)[2],length(NES_TAU$pathway)),
              rep(colnames(NES_TAU)[3],length(NES_TAU$pathway))),
              Model = rep('TAU',length(NES_TAU$pathway)*2),
              Condition =rep(colnames(NES_TAU)[2],length(NES_TAU$pathway)*2),
              Area = NES_TAU$Area,
              padj = c(NES_TAU[,5],NES_TAU[,6]))



nes_mouse <- rbind(nes_tabs[[1]],nes_tabs[[2]],NES_TAU)
nes_mouse <- nes_mouse[nes_mouse$modules != 'Not.Correlated',]

nes_mouse$modules <- factor(nes_mouse$modules,
                            levels = c('M1','M2','M3','M4','M5','M6','M7'))


custom_pal <- c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                "#D6604D", "#B2182B", "#67001F")
custom_pal <- colorRampPalette(custom_pal)(100)



g<-ggplot(nes_mouse, aes_(x=~Group, y=~modules, size=~abs(NES), fill=~NES)) +
  geom_point(color = "white", shape=21) +
  scale_fill_gradientn(colours=custom_pal, space = "Lab") +
  scale_size(range=c(0,30), limits=c(0, NA)) +
  guides(size="none") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 15)) +
  scale_x_discrete(position = "top")



CCX<- g %+% subset(nes_mouse[nes_mouse$padj < 0.01,],Model %in% 'FAD' & Area %in% 'CCX')
HIP<- g %+% subset(nes_mouse[nes_mouse$padj < 0.01,],Model %in% 'FAD' & Area %in% 'HIP')
TAU<- g %+% subset(nes_mouse[nes_mouse$padj < 0.01,],Model %in% 'TAU')

ggsave(file="NEW_RESULTS/results/CCX_NES.jpeg", plot = CCX,units = 'cm' ,width = 15 ,height = 35,
       dpi =1200)

ggsave(file="NEW_RESULTS/results/HIP_NES.jpeg", plot = HIP,units = 'cm' ,width = 15 ,height = 20,
       dpi =1200)

ggsave(file="NEW_RESULTS/results/TAU_HIP_NES.jpeg", plot = TAU,units = 'cm' ,width = 15 ,height = 35,
       dpi =1200)

