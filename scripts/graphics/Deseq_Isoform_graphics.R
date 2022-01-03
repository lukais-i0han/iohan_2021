#### Month and Region analysis
library('dplyr')
library('stringr')
library('ggplot2')
library('cowplot')
library("UpSetR")

## Human

mayo_TCX_DEG <- read.csv('NEW_RESULTS/Deseq_results/human/TCX/mayo_TCX_all.csv',header = T, 
                 stringsAsFactors = F,
                 row.names = 1)

mayo_TCX_DEG <- mayo_TCX_DEG %>% mutate(DEG = case_when(
  DEG == 'DEG' ~ 'DEG',
  TRUE ~ 'Non-Significant'
))

mayo_TCX_DEG <- mayo_TCX_DEG[mayo_TCX_DEG$DEG == 'DEG',]
mayo_TCX_DEG <- mayo_TCX_DEG[mayo_TCX_DEG$Group != 'PA',]   

mayo_TCX_DEG$Group <- ifelse(mayo_TCX_DEG$Group == 'ALZ','AD','PSP')

mayo_TCX_DEG <- mayo_TCX_DEG  %>%
  mutate(UP_DOWN = case_when(
    log2FoldChange <= -log2(1.3) ~ 'Down',
    TRUE ~ 'Up'
  ))

mayo_TCX_DEG <- mayo_TCX_DEG %>% select(Gene_Symbol,log2FoldChange,padj,DEG,
                                  UP_DOWN,Age_Group,Group)

mayo_TCX_DTU <- readRDS('NEW_RESULTS/isoform_results/human/mayo_ISW_TCX_all.rds')

mayo_TCX_DTU <- mayo_TCX_DTU %>% mutate(DTU = case_when(
  DTU == 'DTU' ~ 'gDTU',
  TRUE ~ 'Non-Significant'
))

mayo_TCX_DTU <- mayo_TCX_DTU[mayo_TCX_DTU$DTU == 'gDTU',]
mayo_TCX_DTU <- mayo_TCX_DTU[mayo_TCX_DTU$condition_2 != 'Pathologic_Aging',]
mayo_TCX_DTU <- mayo_TCX_DTU %>% select(gene_name,gene_id,isoform_id,iso_biotype,
                                        dIF,isoform_switch_q_value,Age_Group,Group,DTU)

mayo_TCX_DTU$Group <- ifelse(mayo_TCX_DTU$Group == 'ALZ','AD','PSP')
mayo_TCX_DTU <- mayo_TCX_DTU %>% group_by(Group,Age_Group) %>%
  distinct(gene_name,.keep_all = T)

mayo_DEG_DTU <- list(AD_DEG = mayo_TCX_DEG[mayo_TCX_DEG$Group == 'AD',8],
                     AD_gDTU = mayo_TCX_DTU[mayo_TCX_DTU$Group == 'AD',1],
                     PSP_DEG = mayo_TCX_DEG[mayo_TCX_DEG$Group == 'PSP',8],
                     PSP_gDTU = mayo_TCX_DTU[mayo_TCX_DTU$Group == 'PSP',1])



pdf(file="NEW_RESULTS/results/upset_TCX.pdf",width = 10 ,onefile=FALSE)

upset(fromList(mayo_DEG_DTU),order.by = 'freq',text.scale = 2,line.size = 1,
      group.by = 'degree')

dev.off()

### GGplot2 theme

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

## DEG

j = ggplot2::ggplot(mayo_TCX_DEG, ggplot2::aes(log2FoldChange,-log10(padj))) +
  ggplot2::geom_point(ggplot2::aes(color = DEG),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('red','black')) 

g = j + ggrepel::geom_text_repel(data=mayo_TCX_DEG,point.padding = 20)+
  ggplot2::aes(label="")+
  ggtitle("")  + theme_classic2()+
  geom_vline(xintercept =c(-log2(1.3),log2(1.3)),color = "red",linetype='dashed')+
  geom_hline(yintercept = -log10(0.01),linetype='dashed',color='red')+ 
  theme(
    axis.text.x = element_text(face = "bold", color = "black",size = 20),
    axis.title = element_text(face='bold',color='black',size=20),
    axis.text.y = element_text(face = 'bold',color = 'black',size = 20)
    
  )  +
  facet_wrap(~Age_Group,ncol = 3) + ylim(0,20)+ xlim(c(-8,8))+
  labs(x= expression(log[2](FC)),y= expression(-log[10](FDR))) +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black", face = "bold"
    )
  )


ALZ_DEG<-g %+% subset(mayo_TCX_DEG,Group %in% 'ALZ')
PSP_DEG<-g %+% subset(mayo_TCX_DEG,Group %in% 'PSP')
PA_DEG<-g %+% subset(mayo_TCX_DEG,Group %in% 'PA')


## DTU

h = ggplot2::ggplot(mayo_TCX_DTU, ggplot2::aes(dIF,-log10(isoform_switch_q_value))) +
  ggplot2::geom_point(ggplot2::aes(color = DTU),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('blue','black')) 

i = h + ggrepel::geom_text_repel(data=mayo_TCX_DTU,point.padding = 20)+
  ggplot2::aes(label="")+
  ggtitle("")  + theme_classic2()+ 
  geom_vline(xintercept =c(-0.05,0.05),color = "red",linetype='dashed')+
  geom_hline(yintercept = -log10(0.01),linetype='dashed',color='red')+ 
  theme(
    axis.text.x = element_text(face = "bold", color = "black",size = 20),
    axis.title = element_text(face='bold',color='black',size=20),
    axis.text.y = element_text(face = 'bold',color = 'black',size = 20)
    
  )  +
  facet_wrap(~Age_Group,ncol = 3) + ylim(0,20)+ xlim(c(-0.35,0.35))+
  labs(x= expression(dIF),y= expression(-log[10](FDR))) +
  theme(
    strip.text.x = element_text(
      size = 20, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 20, color = "black", face = "bold"
    )
    )

ALZ_DTU<- i %+% subset(mayo_TCX_DTU,Group %in% 'ALZ')
PSP_DTU<- i %+% subset(mayo_TCX_DTU,Group %in% 'PSP')
PA_DTU<- i %+% subset(mayo_TCX_DTU,Group %in% 'PA')



### PLOTS DEG-DTU per group and age

ALZ <-plot_grid(ALZ_DEG,ALZ_DTU,labels = c('DEG','DTU'),nrow = 2,ncol = 1,label_size = 14)
PSP <-plot_grid(PSP_DEG,PSP_DTU,labels = c('DEG','DTU'),nrow = 2,label_size = 10)
PA <- plot_grid(PA_DEG,PA_DTU,labels = c('DEG','DTU'),nrow = 2,label_size = 10)

ggsave(file="NEW_RESULTS/results/figure2/ALZ.jpeg",plot = ALZ, units = 'cm' ,
       width = 25 ,height = 20,dpi =600)

ggsave(file="NEW_RESULTS/results/figure2/PSP.jpeg",plot = PSP, units = 'cm' ,
       width = 25 ,height = 20,dpi =600)

ggsave(file="NEW_RESULTS/results/figure2/PA.jpeg",plot = PA, units = 'cm' ,
       width = 25 ,height = 20,dpi =600)

### Mouse

FAD5X <- read.csv('NEW_RESULTS/Deseq_results/mouse_FAD5X/FAD5X_all.csv',header = T, 
         stringsAsFactors = F,
         row.names = 1)

FAD5X <- FAD5X %>% mutate(DEG=case_when(
  DEG == 'DEG' ~ 'DEG',
  TRUE ~ 'Non-Significant'
))

FAD5X_HIP <- FAD5X[FAD5X$Area == 'Hippocampus',]
FAD5X_HIP <- FAD5X_HIP[FAD5X_HIP$DEG == 'DEG',]

FAD5X_HIP <- FAD5X_HIP  %>%
  mutate(UP_DOWN = case_when(
    log2FoldChange <= -log2(1.3) ~ 'Down',
    TRUE ~ 'Up'
  ))

FAD5X_HIP <- FAD5X_HIP %>% select(Gene_Symbol,log2FoldChange,padj,DEG,
                                        UP_DOWN,Month)
FAD5X_HIP$'Group' <- '5XFAD'
colnames(FAD5X_HIP) <- c('Gene_Symbol','log2FoldChange','padj','DEG',
                         'UP_DOWN','Age_Group','Group')

FAD5X_DTU <- readRDS('NEW_RESULTS/isoform_results/mouse_FAD5X/FAD5X_ISW.rds')
FAD5X_DTU <- FAD5X_DTU %>% mutate(DTU = case_when(
  DTU == 'DTU' ~ 'gDTU',
  TRUE ~ 'Non-Significant'
))



FAD5X_DTU <- FAD5X_DTU[FAD5X_DTU$Region == 'Hippocampus' &
                       FAD5X_DTU$DTU == 'gDTU',]

FAD5X_DTU <- FAD5X_DTU %>% select(gene_name,gene_id,isoform_id,iso_biotype,
                                  dIF,isoform_switch_q_value,Age_Group,DTU)
FAD5X_DTU$'Group' <- 'FAD5X'

FAD5X_DTU <- FAD5X_DTU %>% group_by(Age_Group) %>%
  distinct(gene_name,.keep_all = T)

Tau35 <- read.csv('NEW_RESULTS/Deseq_results/mouse_TAU/Tau_all.csv',header = T,
         stringsAsFactors = F,
         row.names = 1)

Tau35 <- Tau35 %>% mutate(DEG=case_when(
  DEG == 'DEG' ~ 'DEG',
  TRUE ~ 'Non-Significant'
))

Tau35 <- Tau35[Tau35$DEG == 'DEG',]

Tau35 <- Tau35  %>%
  mutate(UP_DOWN = case_when(
    log2FoldChange <= -log2(1.3) ~ 'Down',
    TRUE ~ 'Up'
  ))

Tau35 <- Tau35 %>% select(Gene_Name,log2FoldChange,padj,DEG,
                                  UP_DOWN,Month)
colnames(Tau35) <- c('Gene_Symbol','log2FoldChange','padj',"DEG",
                     'UP_DOWN','Age_Group')
Tau35$'Group' <- 'TauD35'


Tau <- readRDS('NEW_RESULTS/isoform_results/mouse_TAU/Tau_ISW.rds')
Tau$Age_Group <- factor(Tau$Age_Group,levels = c('4','17'))
Tau <- Tau %>% mutate(DTU = case_when(
  DTU == 'DTU' ~ 'gDTU',
  TRUE ~ 'Non-Significant'
))

Tau <- Tau[Tau$DTU == 'gDTU',]

Tau <- Tau %>% select(gene_name,gene_id,isoform_id,iso_biotype,
                                  dIF,isoform_switch_q_value,Age_Group,DTU)
Tau$'Group' <- 'TauD35'

Tau <- Tau %>% group_by(Age_Group) %>%
  distinct(gene_name,.keep_all = T)

fad5x_tau <- list("5XFAD_DEG" = FAD5X_HIP$Gene_Symbol,
                  "5XFAD_gDTU" = FAD5X_DTU$gene_name,
                  "TauD35_DEG" = Tau35$Gene_Symbol,
                  "TauD35_gDTU" = Tau$gene_name)


pdf(file="NEW_RESULTS/results/upset_mouse.pdf",width = 10 ,onefile=FALSE)

upset(fromList(fad5x_tau),order.by = 'freq',text.scale = 2,line.size = 1,
      group.by = 'degree')

dev.off()

human_mouse_DEG <- rbind(mayo_TCX_DEG,FAD5X_HIP,Tau35)
write.csv(human_mouse_DEG,'~/Desktop/Mestrado/Paper-Lukas/tables/human_mouse_DEG.csv')

human_mouse_DTU <- rbind(mayo_TCX_DTU,FAD5X_DTU,Tau)
write.csv(human_mouse_DTU,'~/Desktop/Mestrado/Paper-Lukas/tables/human_mouse_DTU.csv')

#### FAD-DEG

FAD5X$Month <- factor(FAD5X$Month, levels = c('4','12','18'))


j = ggplot2::ggplot(FAD5X, ggplot2::aes(log2FoldChange,-log10(padj))) +
  ggplot2::geom_point(ggplot2::aes(color = DEG),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('red','black')) 

g = j + ggrepel::geom_text_repel(data=FAD5X,point.padding = 20)+
  ggplot2::aes(label="")+
  ggtitle("")  + theme_classic2()+
  geom_vline(xintercept =c(-log2(1.3),log2(1.3)),color = "red",linetype='dashed')+
  geom_hline(yintercept = -log10(0.01),linetype='dashed',color='red')+ 
  theme(
        axis.text.x = element_text(face = "bold", color = "black",size = 15),
        axis.title = element_text(face='bold',color='black',size=15),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
        
  )  +
  facet_wrap(~Month,ncol = 3) + ylim(0,20)+ xlim(c(-4.5,4.5))+
  labs(x= expression(log[2](FC)),y= expression(-log[10](FDR))) +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 15, color = "black", face = "bold"
    )
  )


CCX<-g %+% subset(FAD5X,Area %in% 'Cortex')
HIP<-g %+% subset(FAD5X,Area %in% 'Hippocampus')


FAD<-plot_grid(CCX,HIP,labels = c('Cortex','Hippocampus'),nrow = 2,label_size = 10)

ggsave(file="NEW_RESULTS/results/figure4/FAD_DEG.jpeg",plot = FAD, units = 'cm' ,width = 30 ,height = 15,
       dpi =600)


##DTU

FAD5X_DTU$Age_Group <- factor(FAD5X_DTU$Age_Group, levels = c('4','12','18'))

h = ggplot2::ggplot(FAD5X_DTU, ggplot2::aes(dIF,-log10(iso_q_value))) +
  ggplot2::geom_point(ggplot2::aes(color = DTU),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('blue','black')) 

i = h + ggrepel::geom_text_repel(data=FAD5X_DTU,point.padding = 20)+
  ggplot2::aes(label="")+
  ggtitle("")  + theme_classic2()+
  geom_vline(xintercept =c(-0.05,0.05),color = "red",linetype='dashed')+
  geom_hline(yintercept = -log10(0.01),linetype='dashed',color='red')+ 
  theme(legend.position = 'right',
        axis.text.x = element_text(face = "bold", color = "black",size = 15),
        axis.title = element_text(face='bold',color='black',size=15),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
        
  )  +
  theme(
        legend.box.just = 'right')+
  facet_wrap(~Age_Group,ncol = 3) + ylim(0,20)+ xlim(c(-0.35,0.35))+
  labs(x= expression(dIF),y= expression(-log[10](FDR))) +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 15, color = "black", face = "bold"
    )
  )

CCX_DTU<- i %+% subset(FAD5X_DTU,Region %in% 'Cortex')
HIP_DTU<-i %+% subset(FAD5X_DTU,Region %in% 'Hippocampus')

FAD_CCX<- plot_grid(CCX,CCX_DTU,labels = c('DEG','gDTU'),nrow = 2,label_size = 10)
FAD_HIP<- plot_grid(HIP,HIP_DTU,labels = c('DEG','gDTU'),nrow = 2,label_size = 10)

ggsave(file="NEW_RESULTS/results/figure4/FAD_CCX.jpeg",plot = FAD_CCX, units = 'cm' ,width = 30 ,height = 15,
       dpi =600)

ggsave(file="NEW_RESULTS/results/figure4/FAD_HIP.jpeg",plot = FAD_HIP, units = 'cm' ,width = 30 ,height = 15,
       dpi =600)
### TAU

##DEG
Tau35$Month <- factor(Tau35$Month,levels = c('4','17'))

j = ggplot2::ggplot(Tau35, ggplot2::aes(log2FoldChange,-log10(padj))) +
  ggplot2::geom_point(ggplot2::aes(color = DEG),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('red','black')) 

g = j + ggrepel::geom_text_repel(data=,point.padding = 20)+
  ggplot2::aes(label="")+
  ggtitle("")  + theme_classic2()+
  geom_vline(xintercept =c(-log2(1.3),log2(1.3)),color = "red",linetype='dashed')+
  geom_hline(yintercept = -log10(0.01),linetype='dashed',color='red')+ 
  theme(
        axis.text.x = element_text(face = "bold", color = "black",size = 15),
        axis.title = element_text(face='bold',color='black',size=15),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
        
  )  +
  theme(
        legend.box.just = 'right')+
  facet_wrap(~Month,ncol = 2) + ylim(0,20)+ xlim(c(-8,8))+
  labs(x= expression(log[2](FC)),y= expression(-log[10](padj))) +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 15, color = "black", face = "bold"
    )
  )
g
## DTU

Tau$Age_Group <- factor(Tau$Age_Group,levels = c('4','17'))

h = ggplot2::ggplot(Tau, ggplot2::aes(dIF,-log10(iso_q_value))) +
  ggplot2::geom_point(ggplot2::aes(color = DTU),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('blue','black')) 

i = h + ggrepel::geom_text_repel(data=Tau,point.padding = 20)+
  ggplot2::aes(label="")+
  ggtitle("")  + theme_classic2()+
  geom_vline(xintercept =c(-0.05,0.05),color = "red",linetype='dashed')+
  geom_hline(yintercept = -log10(0.01),linetype='dashed',color='red')+ 
  theme(legend.position = 'right',
        axis.text.x = element_text(face = "bold", color = "black",size = 15),
        axis.title = element_text(face='bold',color='black',size=15),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
        
  )  +
  theme(legend.justification = 'top',
        legend.box.just = 'right')+
  facet_wrap(~Age_Group,ncol = 2) + ylim(0,20)+ xlim(c(-0.35,0.35))+
  labs(x= expression(dIF),y= expression(-log[10](FDR))) +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 15, color = "black", face = "bold"
    )
  )

TAU<-plot_grid(g,i,labels = c('DEG','DTU'),nrow = 2,label_size = 10)

ggsave(file="NEW_RESULTS/results/figure5/TAU.jpeg",plot = TAU, 
       units = 'cm' ,width = 20 ,height = 25,
       dpi =600)


