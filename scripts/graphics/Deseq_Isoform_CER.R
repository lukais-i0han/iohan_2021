library('dplyr')
library('stringr')
library('ggplot2')
library('cowplot')

## Human

## Human

mayo_CER_DEG <- read.csv('NEW_RESULTS/Deseq_results/human/CER/mayo_CER_all.csv',header = T, 
                         stringsAsFactors = F,
                         row.names = 1)

mayo_CER_DEG <- mayo_CER_DEG %>% mutate(DEG = case_when(
  DEG == 'DEG' ~ 'DEG',
  TRUE ~ 'Non-Significant'
))

mayo_CER_DEG <- mayo_CER_DEG[mayo_CER_DEG$DEG == 'DEG',]
mayo_CER_DEG <- mayo_CER_DEG[mayo_CER_DEG$Group != 'PA',]   

mayo_CER_DEG$Group <- ifelse(mayo_CER_DEG$Group == 'ALZ','AD','PSP')

mayo_CER_DEG <- mayo_CER_DEG  %>%
  mutate(UP_DOWN = case_when(
    log2FoldChange <= -log2(1.3) ~ 'Down',
    TRUE ~ 'Up'
  ))

mayo_CER_DEG <- mayo_CER_DEG %>% select(Gene_Symbol,log2FoldChange,padj,DEG,
                                        UP_DOWN,Age_Group,Group)

mayo_CER_DTU <- readRDS('NEW_RESULTS/isoform_results/human/mayo_ISW_CER_all.rds')

mayo_CER_DTU <- mayo_CER_DTU %>% mutate(DTU = case_when(
  DTU == 'DTU' ~ 'gDTU',
  TRUE ~ 'Non-Significant'
))

mayo_CER_DTU <- mayo_CER_DTU[mayo_CER_DTU$DTU == 'gDTU',]
mayo_CER_DTU <- mayo_CER_DTU[mayo_CER_DTU$condition_2 != 'Pathologic_Aging',]
mayo_CER_DTU <- mayo_CER_DTU %>% select(gene_name,gene_id,isoform_id,iso_biotype,
                                        dIF,isoform_switch_q_value,Age_Group,Group,DTU)

mayo_CER_DTU$Group <- ifelse(mayo_CER_DTU$Group == 'ALZ','AD','PSP')

write.csv(mayo_CER_DEG,'~/Desktop/Mestrado/Paper-Lukas/tables/human_CER_DEG.csv')
write.csv(mayo_CER_DTU,'~/Desktop/Mestrado/Paper-Lukas/tables/human_CER_DTU.csv')

mayo_DEG_DTU <- list(AD_DEG = mayo_CER_DEG[mayo_CER_DEG$Group == 'AD',1],
                     AD_gDTU = mayo_CER_DTU[mayo_CER_DTU$Group == 'AD',1],
                     PSP_DEG = mayo_CER_DEG[mayo_CER_DEG$Group == 'PSP',1],
                     PSP_gDTU = mayo_CER_DTU[mayo_CER_DTU$Group == 'PSP',1])




pdf(file="NEW_RESULTS/results/upset_CER.pdf",width = 10 ,onefile=FALSE)

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

j = ggplot2::ggplot(mayo_CER_DEG, ggplot2::aes(log2FoldChange,-log10(padj))) +
  ggplot2::geom_point(ggplot2::aes(color = DEG),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('red','black')) 

g = j + ggrepel::geom_text_repel(data=mayo_CER_DEG,point.padding = 20)+
  ggplot2::aes(label="")+
  ggtitle("")  + theme_classic2()+
  geom_vline(xintercept =c(-log2(1.3),log2(1.3)),color = "red",linetype='dashed')+
  geom_hline(yintercept = -log10(0.01),linetype='dashed',color='red')+ 
  theme(
    axis.text.x = element_text(face = "bold", color = "black",size = 15),
    axis.title = element_text(face='bold',color='black',size=15),
    axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
    
  )  +
  facet_wrap(~Age_Group,ncol = 3) + ylim(0,20)+ xlim(c(-8,8))+
  labs(x= expression(log[2](FC)),y= expression(-log[10](FDR))) +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 15, color = "black", face = "bold"
    )
  )


ALZ_DEG<-g %+% subset(mayo_CER_DEG,Group %in% 'ALZ')
PSP_DEG<-g %+% subset(mayo_CER_DEG,Group %in% 'PSP')
PA_DEG<-g %+% subset(mayo_CER_DEG,Group %in% 'PA')


## DTU

h = ggplot2::ggplot(mayo_CER_DTU, ggplot2::aes(dIF,-log10(isoform_switch_q_value))) +
  ggplot2::geom_point(ggplot2::aes(color = DTU),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('blue','black')) 

i = h + ggrepel::geom_text_repel(data=mayo_CER_DTU,point.padding = 20)+
  ggplot2::aes(label="")+
  ggtitle("")  + theme_classic2()+ 
  geom_vline(xintercept =c(-0.05,0.05),color = "red",linetype='dashed')+
  geom_hline(yintercept = -log10(0.01),linetype='dashed',color='red')+ 
  theme(
    axis.text.x = element_text(face = "bold", color = "black",size = 15),
    axis.title = element_text(face='bold',color='black',size=15),
    axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
    
  )  +
  facet_wrap(~Age_Group,ncol = 3) + ylim(0,20)+ xlim(c(-0.35,0.35))+
  labs(x= expression(dIF),y= '') +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_blank()
  )

ALZ_DTU<- i %+% subset(mayo_CER_DTU,Group %in% 'ALZ')
PSP_DTU<- i %+% subset(mayo_CER_DTU,Group %in% 'PSP')
PA_DTU<- i %+% subset(mayo_CER_DTU,Group %in% 'PA')



### PLOTS DEG-DTU per group and age

ALZ <-plot_grid(ALZ_DEG,ALZ_DTU,labels = c('DEG','DTU'),nrow = 2,label_size = 10)
PSP <-plot_grid(PSP_DEG,PSP_DTU,labels = c('DEG','DTU'),nrow = 2,label_size = 10)
PA <- plot_grid(PA_DEG,PA_DTU,labels = c('DEG','DTU'),nrow = 2,label_size = 10)

ggsave(file="NEW_RESULTS/results/sup1//ALZ_CER.jpeg",plot = ALZ, units = 'cm' ,
       width = 25 ,height = 20,dpi =600)

ggsave(file="NEW_RESULTS/results/figure2/PSP_CER.jpeg",plot = PSP, units = 'cm' ,
       width = 25 ,height = 20,dpi =600)

ggsave(file="NEW_RESULTS/results/figure2/PA_CER.jpeg",plot = PA, units = 'cm' ,
       width = 25 ,height = 20,dpi =600)
