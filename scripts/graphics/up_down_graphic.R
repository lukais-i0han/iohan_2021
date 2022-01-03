### Plot of groups
library('ggplot2')
library('dplyr')

## Mayo TCX

mayo_TCX <- read.csv('NEW_RESULTS/Deseq_results/human/TCX/mayo_TCX_sig.csv',header = T,
                     stringsAsFactors = T, row.names = 1)


mayo_TCX <- mayo_TCX  %>%
  mutate(DEG = case_when(
    log2FoldChange <= -log2(1.3) ~ 'Down',
    TRUE ~ 'Up'
  ))


mayo_plot <- ggplot(mayo_TCX,aes(x=Age_Group,
          y=factor(Gene_Symbol, levels = unique((Gene_Symbol)))))+
  geom_jitter(aes(color=DEG),size=2)+
  labs(y='Genes', x=NULL,fill = 'DEG') + 
  scale_color_manual(values = c('firebrick','blue'),labels = c('Up','Down'))+
  theme_classic2(base_size = 10) + theme(strip.background = element_blank(),
                                     axis.text.x = element_text(face = "bold", color = "black",size = 12),
                                     axis.text.y = element_blank(),
                                     legend.text = element_text(size = 12,face = 'bold'),
                                     legend.justification = 'top',
                                     legend.box.just = 'right') +facet_grid(~Group)


ALZ_TCX <- mayo_plot %+% subset(mayo_TCX, Group %in% 'ALZ') 
PSP_TCX <- mayo_plot %+% subset(mayo_TCX, Group %in% 'PSP') 
PA_TCX <- mayo_plot %+% subset(mayo_TCX, Group %in% 'PA') 


ggsave('NEW_RESULTS/results/ALZ_TCX.jpeg',plot = ALZ_TCX,height = 15,
       width = 20,units = 'cm',dpi=1200)
ggsave('NEW_RESULTS/results/PSP_TCX.jpeg',plot = PSP_TCX,height = 15,
       width = 20,units = 'cm',dpi=1200)
ggsave('NEW_RESULTS/results/PA_TCX.jpeg',plot = PA_TCX,height = 15,
       width = 20,units = 'cm',dpi=1200)

#### CER

mayo_CER <- read.csv('NEW_RESULTS/Deseq_results/human/CER/mayo_CER_sig.csv',header = T, 
                    stringsAsFactors = F,row.names = 1)


mayo_CER<- mayo_CER %>% select(log2FoldChange,Gene_Symbol,DEG,Age_Group,Group) %>%
  mutate(DEG = case_when(
    log2FoldChange <= -log2(1.3) ~ 'Down',
    TRUE ~ 'Up'
  ))

mayo_plot_CER <- ggplot(mayo_CER,aes(x=Age_Group,
                                 y=factor(Gene_Symbol, levels = unique((Gene_Symbol)))))+
  geom_jitter(aes(color=DEG),size=2)+
  labs(y='Genes', x=NULL,fill = 'DEG') + 
  scale_color_manual(values = c('firebrick','blue'),labels = c('Up','Down'))+
  theme_classic2(base_size = 10) + theme(strip.background = element_blank(),
                                         axis.text.x = element_text(face = "bold", color = "black",size = 12),
                                         axis.text.y = element_blank(),
                                         legend.text = element_text(size = 12,face = 'bold'),
                                         legend.justification = 'top',
                                         legend.box.just = 'right') +facet_grid(~Group)



ALZ_CER <- mayo_plot_CER %+% subset(mayo_CER, Group %in% 'ALZ') 
PSP_CER <- mayo_plot_CER %+% subset(mayo_CER, Group %in% 'PSP') 
PA_CER <- mayo_plot_CER %+% subset(mayo_CER, Group %in% 'PA') 


ggsave('NEW_RESULTS/results/ALZ_CER.jpeg',plot = ALZ_CER,height = 15,
      width = 20,units = 'cm',dpi=1200)
ggsave('NEW_RESULTS/results/PSP_CER.jpeg',plot = PSP_CER,height = 15,
       width = 20,units = 'cm',dpi=1200)
ggsave('NEW_RESULTS/results/PA_CER.jpeg',plot = PA_CER,height = 15,
       width = 20,units = 'cm',dpi=1200)
#### FAD-5X

FAD5X <- read.csv('Deseq_results/mouse_FAD5X/FAD5X_all_sig.csv',row.names = 1,
                  stringsAsFactors = F)

FAD5X <- FAD5X %>% select(log2FoldChange,Gene_Symbol,Month,Area,DEG) %>%
  mutate(DEG = case_when(
      log2FoldChange <= -log2(1.3) ~ 'Down',
      TRUE ~ 'Up'
  ))

FAD5X <- FAD5X %>% mutate(Month = case_when(
  Month == '4' ~ '[4]',
  Month == '12' ~ '[12]',
  TRUE ~ '[18]'
))

FAD5X$Month <- factor(FAD5X$Month,levels = c('[4]','[12]','[18]'))
FAD5X <- FAD5X[FAD5X$Area == 'Hippocampus',]

myplot_FAD <-ggplot(FAD5X,aes(x=Month,
                                 y=factor(Gene_Symbol, levels = unique((Gene_Symbol)))))+
  geom_jitter(aes(color=DEG),size=2)+ 
  labs(y='Genes', x=NULL,fill = 'DEG') + 
  scale_color_manual(values = c('firebrick','blue'),labels = c('Up','Down'))+
  theme_classic2(base_size = 10) + theme(strip.background = element_blank(),
                                     axis.text.x = element_text(face = "bold", color = "black",size = 12),
                                     axis.text.y = element_blank(),
                                     legend.text = element_text(size = 12,face = 'bold'),
                                     legend.justification = 'top',
                                     legend.box.just = 'right')  

ggsave('NEW_RESULTS/results/FAD5X_HIP.jpeg',plot = myplot_FAD,height = 15,width = 10,units = 'cm',dpi=1200)
### TauD35

Taud35 <- read.csv('Deseq_results/mouse_TAU/Tau_sig.csv',row.names = 1,
                   stringsAsFactors = F)

Taud35 <- Taud35 %>% select(log2FoldChange,Gene_Name,Month,DEG) %>%
  mutate(DEG = case_when(
    log2FoldChange <=  -log2(1.3) ~ 'Down',
    TRUE ~ 'Up'
  )) %>% mutate(Month = case_when(
    Month == '4' ~ '[4]',
    TRUE ~ '[17]'
  ))

Taud35$Month <- factor(Taud35$Month,levels = c('[4]','[17]'))

myplot_TAU <-ggplot(Taud35,aes(x=Month,
                              y=factor(Gene_Name, levels = unique((Gene_Name)))))+
  geom_jitter(aes(color=DEG),size=2)+
  labs(y='Genes', x=NULL,fill = 'DEG') + 
  scale_color_manual(values = c('firebrick','blue'),labels = c('Up','Down'))+
  theme_classic2(base_size = 10) + theme(strip.background = element_blank(),
                                     axis.text.x = element_text(face = "bold", color = "black",size = 12),
                                     axis.text.y = element_blank(),
                                     legend.text = element_text(size = 12,face = 'bold'),
                                     legend.justification = 'top',
                                     legend.box.just = 'right')


plots <- list(mayo_plot,myplot_FAD)
plot_name <- c('Human_TCX','FAD5X')
for(i in 1:length(plots)){
  ggsave(file=paste0('NEW_RESULTS/results/plot_',plot_name[i],'.jpeg'),
         plot = plots[[i]],height = 15,width = 25,units = 'cm',dpi = 1200)
}

ggsave('NEW_RESULTS/results/Tau.jpeg',plot = myplot_TAU,height = 15,width = 10,units = 'cm',dpi=1200)


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
