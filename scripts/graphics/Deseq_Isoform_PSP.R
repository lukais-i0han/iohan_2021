mayo_TCX_DEG <- read.csv('NEW_RESULTS/Deseq_results/human/TCX/mayo_TCX_all.csv',header = T, 
                         stringsAsFactors = F,
                         row.names = 1)

mayo_TCX_DEG <- mayo_TCX_DEG %>% mutate(DEG = case_when(
  DEG == 'DEG' ~ 'DEG',
  TRUE ~ 'Non-Significant'
))
mayo_TCX_DEG <- mayo_TCX_DEG[mayo_TCX_DEG$Group == 'PSP',]

mayo_TCX_DTU <- readRDS('NEW_RESULTS/isoform_results/human/mayo_ISW_TCX_all.rds')

mayo_TCX_DTU <- mayo_TCX_DTU %>% mutate(DTU = case_when(
  DTU == 'DTU' ~ 'gDTU',
  TRUE ~ 'Non-Significant'
))
mayo_TCX_DTU <- mayo_TCX_DTU[mayo_TCX_DTU$Group == 'PSP',]



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
PSP <-plot_grid(g,i,labels = c('DEG','DTU'),nrow = 2,label_size = 10)

ggsave(file="NEW_RESULTS/results/figure2/PSP.jpeg",plot = PSP, units = 'cm' ,
       width = 25 ,height = 20,dpi =600)
