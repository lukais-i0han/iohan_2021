library('stringr')
library('readr')
library('RedeR')
library('RColorBrewer')
library('igraph')
library('dplyr')
library('CEMiTool')
library('ggplot2')
library('ggrepel')
library('gprofiler2')

#### PPI cemitool


plot_interaction <- function(ig_obj, n, color, name, mod_genes, coexp_hubs,AD_rik){
  degrees <- igraph::degree(ig_obj, normalized=FALSE)
  ig_obj <- igraph::set_vertex_attr(ig_obj, "degree", value = degrees)
  max_n <- min(n, length(degrees))
  net_obj <- intergraph::asNetwork(ig_obj)
  m <- network::as.matrix.network.adjacency(net_obj) # get sociomatrix
  # get coordinates from Fruchterman and Reingold's force-directed placement algorithm.
  plotcord <- data.frame(sna::gplot.layout.fruchtermanreingold(m, NULL))
  # or get it them from Kamada-Kawai's algorithm:
  # plotcord <- data.frame(sna::gplot.layout.kamadakawai(m, NULL))
  colnames(plotcord) <- c("X1","X2")
  edglist <- network::as.matrix.network.edgelist(net_obj)
  edges <- data.frame(plotcord[edglist[,1],], plotcord[edglist[,2],])
  plotcord$vertex.names <- as.factor(network::get.vertex.attribute(net_obj, "vertex.names"))
  plotcord$Degree <- network::get.vertex.attribute(net_obj, "degree")
  plotcord[, "shouldLabel"] <- FALSE
  plotcord[, "Hub"] <- ""
  int_hubs <- names(sort(degrees, decreasing=TRUE))[1:max_n]
  int_bool <- plotcord[, "vertex.names"] %in% int_hubs
  plotcord[which(int_bool), "Hub"] <- "Interaction"
  sel_vertex <- int_hubs
  if(!missing(coexp_hubs)){
    coexp_bool <- plotcord[, "vertex.names"] %in% coexp_hubs
    coexp_and_int <- coexp_bool & int_bool
    plotcord[which(coexp_bool), "Hub"] <- "Co-expression"
    plotcord[which(coexp_and_int), "Hub"] <- "Co-expression + Interaction"
    sel_vertex <- c(sel_vertex, coexp_hubs)
  }
  
  colnames(edges) <-  c("X1","Y1","X2","Y2")
  #edges$midX  <- (edges$X1 + edges$X2) / 2
  #edges$midY  <- (edges$Y1 + edges$Y2) / 2
  plotcord[which(plotcord[, "vertex.names"] %in% sel_vertex), "shouldLabel"] <- TRUE
  plotcord$Degree_cut <- cut(plotcord$Degree, breaks=3, labels=FALSE)
  plotcord$in_mod <- TRUE
  #mod_genes <- cem@module[cem@module$modules==name,]$genes
  not_in <- setdiff(plotcord[,"vertex.names"], mod_genes)
  plotcord[which(plotcord[, "vertex.names"] %in% not_in), "in_mod"] <- FALSE
  plotcord$Hub <- ifelse(plotcord$vertex.names %in% AD_risk$name,'AD risk gene',plotcord$Hub)
  plotcord$shouldLabel <- ifelse(plotcord$Hub == 'AD risk gene' &
                                   plotcord$in_mod == TRUE,TRUE,plotcord$shouldLabel)
  
  pl <- ggplot(plotcord)  +
    geom_segment(data=edges, aes_(x=~X1, y=~Y1, xend=~X2, yend=~Y2),
                 size = 0.5, alpha=0.5, colour="#DDDDDD") +
    geom_point(aes_(x=~X1, y=~X2, size=~Degree, alpha=~Degree), color=color) +
    geom_label_repel(aes_(x=~X1, y=~X2, label=~vertex.names, color=~Hub),
                     box.padding=unit(1, "lines"),
                     data=function(x){x[x$shouldLabel, ]}) +
    scale_colour_manual(values=c("Co-expression" = "#005E87",
                                 "Interaction" = "#540814",
                                 "Co-expression + Interaction" = "#736E0B",
                                 "AD risk gene" = 'black')) +
    labs(title=name) +
    ggplot2::theme_bw(base_size = 12, base_family = "") +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   legend.key = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(fill = "white",
                                                            colour = NA),
                   panel.border = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
  
  return(pl)
}


human_cem_ALZ <- readRDS('NEW_RESULTS/Cemtool_results/human/pval_01/ALZ.rds')

AD_risk <- read.csv('NEW_RESULTS/ad_risk.csv',header = T, stringsAsFactors = F)
AD_risk <- AD_risk[str_detect(AD_risk$genes,'ENSG'),] 

AD_risk <- gconvert(AD_risk,organism = 'hsapiens',filter_na = T)
AD_risk <- AD_risk[AD_risk$name != 'nan',]

mod_cols <- mod_colors(human_cem_ALZ)
mod_names <- names(human_cem_ALZ@interactions)
mod_names <- mod_names[which(mod_names!="Not.Correlated")]
hubs <- lapply(get_hubs(human_cem_ALZ), names)
zero_ints <- character()
zero_ints <- lapply(names(human_cem_ALZ@interactions), function(mod){
  degree <- igraph::degree(human_cem_ALZ@interactions[[mod]], normalized=FALSE)
  if(length(degree) == 0) {
    zero_ints <- append(zero_ints, mod)
  }
})
zero_ints <- unlist(zero_ints)
if(!is.null(zero_ints)){
  mod_names <- mod_names[which(!(mod_names %in% zero_ints))]
}

modules_colors <- c('#F8766D','#CD9600','#7CAE00',
                    '#00BE67','#00BFC4','#00A9FF',
                    '#C77CFF','#FF61CC')
names(modules_colors) <- mod_names

res <- lapply(mod_names, function(name){
  plot_interaction(ig_obj=human_cem_ALZ@interactions[[name]],
                   n=10, color=modules_colors[name], name=name,
                   mod_genes=module_genes(human_cem_ALZ, name)$genes,
                   coexp_hubs=hubs[[name]],
                   AD_rik = AD_risk)
})
names(res) <- mod_names

plots_ALZ <- c('M1','M2','M3','M4','M5','M6','M7','M8')

for(i in 1:length(plots_ALZ)){
  ggsave(file=paste0('NEW_RESULTS/results/cemitool_01/ppi_plots/ALZ//plot_',
         plots_ALZ[i],'.jpeg'),
         plot = res[[plots_ALZ[i]]],height = 25,width = 30,units = 'cm',dpi =600)
}


### PSP

human_cem_PSP <- readRDS('NEW_RESULTS/Cemtool_results/human/pval_01/PSP.rds')

PSP_risk <- read.csv('NEW_RESULTS/PSP_risk_genes.csv',header = T, stringsAsFactors = F)



mod_cols <- mod_colors(human_cem_PSP)
mod_names <- names(human_cem_PSP@interactions)
mod_names <- mod_names[which(mod_names!="Not.Correlated")]
hubs <- lapply(get_hubs(human_cem_PSP), names)
zero_ints <- character()
zero_ints <- lapply(names(human_cem_PSP@interactions), function(mod){
  degree <- igraph::degree(human_cem_PSP@interactions[[mod]], normalized=FALSE)
  if(length(degree) == 0) {
    zero_ints <- append(zero_ints, mod)
  }
})
zero_ints <- unlist(zero_ints)
if(!is.null(zero_ints)){
  mod_names <- mod_names[which(!(mod_names %in% zero_ints))]
}

plot_interaction <- function(ig_obj, n, color, name, mod_genes, coexp_hubs,PSP_risk){
  degrees <- igraph::degree(ig_obj, normalized=FALSE)
  ig_obj <- igraph::set_vertex_attr(ig_obj, "degree", value = degrees)
  max_n <- min(n, length(degrees))
  net_obj <- intergraph::asNetwork(ig_obj)
  m <- network::as.matrix.network.adjacency(net_obj) # get sociomatrix
  # get coordinates from Fruchterman and Reingold's force-directed placement algorithm.
  plotcord <- data.frame(sna::gplot.layout.fruchtermanreingold(m, NULL))
  # or get it them from Kamada-Kawai's algorithm:
  # plotcord <- data.frame(sna::gplot.layout.kamadakawai(m, NULL))
  colnames(plotcord) <- c("X1","X2")
  edglist <- network::as.matrix.network.edgelist(net_obj)
  edges <- data.frame(plotcord[edglist[,1],], plotcord[edglist[,2],])
  plotcord$vertex.names <- as.factor(network::get.vertex.attribute(net_obj, "vertex.names"))
  plotcord$Degree <- network::get.vertex.attribute(net_obj, "degree")
  plotcord[, "shouldLabel"] <- FALSE
  plotcord[, "Hub"] <- ""
  int_hubs <- names(sort(degrees, decreasing=TRUE))[1:max_n]
  int_bool <- plotcord[, "vertex.names"] %in% int_hubs
  plotcord[which(int_bool), "Hub"] <- "Interaction"
  sel_vertex <- int_hubs
  if(!missing(coexp_hubs)){
    coexp_bool <- plotcord[, "vertex.names"] %in% coexp_hubs
    coexp_and_int <- coexp_bool & int_bool
    plotcord[which(coexp_bool), "Hub"] <- "Co-expression"
    plotcord[which(coexp_and_int), "Hub"] <- "Co-expression + Interaction"
    sel_vertex <- c(sel_vertex, coexp_hubs)
  }
  
  colnames(edges) <-  c("X1","Y1","X2","Y2")
  #edges$midX  <- (edges$X1 + edges$X2) / 2
  #edges$midY  <- (edges$Y1 + edges$Y2) / 2
  plotcord[which(plotcord[, "vertex.names"] %in% sel_vertex), "shouldLabel"] <- TRUE
  plotcord$Degree_cut <- cut(plotcord$Degree, breaks=3, labels=FALSE)
  plotcord$in_mod <- TRUE
  #mod_genes <- cem@module[cem@module$modules==name,]$genes
  not_in <- setdiff(plotcord[,"vertex.names"], mod_genes)
  plotcord[which(plotcord[, "vertex.names"] %in% not_in), "in_mod"] <- FALSE
  plotcord$Hub <- ifelse(plotcord$vertex.names %in% PSP_risk$genes,'PSP risk gene',plotcord$Hub)
  plotcord$shouldLabel <- ifelse(plotcord$Hub == 'PSP risk gene' &
                                   plotcord$in_mod == TRUE,TRUE,plotcord$shouldLabel)
  
  pl <- ggplot(plotcord)  +
    geom_segment(data=edges, aes_(x=~X1, y=~Y1, xend=~X2, yend=~Y2),
                 size = 0.5, alpha=0.5, colour="#DDDDDD") +
    geom_point(aes_(x=~X1, y=~X2, size=~Degree, alpha=~Degree), color=color) +
    geom_label_repel(aes_(x=~X1, y=~X2, label=~vertex.names, color=~Hub),
                     box.padding=unit(1, "lines"),
                     data=function(x){x[x$shouldLabel, ]}) +
    scale_colour_manual(values=c("Co-expression" = "#005E87",
                                 "Interaction" = "#540814",
                                 "Co-expression + Interaction" = "#736E0B",
                                 "PSP risk gene" = 'black')) +
    labs(title=name) +
    ggplot2::theme_bw(base_size = 12, base_family = "") +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   legend.key = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(fill = "white",
                                                            colour = NA),
                   panel.border = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
  
  return(pl)
}

modules_colors <- c('#F8766D','#CD9600','#7CAE00',
                    '#00BE67','#00BFC4','#00A9FF',
                    '#C77CFF','#FF61CC')
names(modules_colors) <- mod_names

res <- lapply(mod_names, function(name){
  plot_interaction(ig_obj=human_cem_PSP@interactions[[name]],
                   n=10, color=modules_colors[name], name=name,
                   mod_genes=module_genes(human_cem_PSP, name)$genes,
                   coexp_hubs=hubs[[name]],
                   PSP_risk = PSP_risk)
})
names(res) <- mod_names

plots_PSP <- c('M1','M2','M3','M4','M5','M6','M7','M8')

for(i in 1:length(plots_PSP)){
  ggsave(file=paste0('NEW_RESULTS/results/cemitool_01/ppi_plots/PSP/plot_',
                     plots_PSP[i],'.jpeg'),
         plot = res[[plots_PSP[i]]],height = 25,width = 30,units = 'cm',dpi = 600)
}


### MSBB 10 and 36

human_cem_MSBB <- readRDS('NEW_RESULTS/Cemtool_results/MSBB/ALZ_10.rds')

AD_risk <- read.csv('NEW_RESULTS/ad_risk.csv',header = T, stringsAsFactors = F)

mod_cols <- mod_colors(human_cem_MSBB)
mod_names <- names(human_cem_MSBB@interactions)
mod_names <- mod_names[which(mod_names!="Not.Correlated")]
hubs <- lapply(get_hubs(human_cem_MSBB), names)
zero_ints <- character()
zero_ints <- lapply(names(human_cem_MSBB@interactions), function(mod){
  degree <- igraph::degree(human_cem_MSBB@interactions[[mod]], normalized=FALSE)
  if(length(degree) == 0) {
    zero_ints <- append(zero_ints, mod)
  }
})
zero_ints <- unlist(zero_ints)
if(!is.null(zero_ints)){
  mod_names <- mod_names[which(!(mod_names %in% zero_ints))]
}

plot_interaction <- function(ig_obj, n, color, name, mod_genes, coexp_hubs,MSBB_risk){
  degrees <- igraph::degree(ig_obj, normalized=FALSE)
  ig_obj <- igraph::set_vertex_attr(ig_obj, "degree", value = degrees)
  max_n <- min(n, length(degrees))
  net_obj <- intergraph::asNetwork(ig_obj)
  m <- network::as.matrix.network.adjacency(net_obj) # get sociomatrix
  # get coordinates from Fruchterman and Reingold's force-directed placement algorithm.
  plotcord <- data.frame(sna::gplot.layout.fruchtermanreingold(m, NULL))
  # or get it them from Kamada-Kawai's algorithm:
  # plotcord <- data.frame(sna::gplot.layout.kamadakawai(m, NULL))
  colnames(plotcord) <- c("X1","X2")
  edglist <- network::as.matrix.network.edgelist(net_obj)
  edges <- data.frame(plotcord[edglist[,1],], plotcord[edglist[,2],])
  plotcord$vertex.names <- as.factor(network::get.vertex.attribute(net_obj, "vertex.names"))
  plotcord$Degree <- network::get.vertex.attribute(net_obj, "degree")
  plotcord[, "shouldLabel"] <- FALSE
  plotcord[, "Hub"] <- ""
  int_hubs <- names(sort(degrees, decreasing=TRUE))[1:max_n]
  int_bool <- plotcord[, "vertex.names"] %in% int_hubs
  plotcord[which(int_bool), "Hub"] <- "Interaction"
  sel_vertex <- int_hubs
  if(!missing(coexp_hubs)){
    coexp_bool <- plotcord[, "vertex.names"] %in% coexp_hubs
    coexp_and_int <- coexp_bool & int_bool
    plotcord[which(coexp_bool), "Hub"] <- "Co-expression"
    plotcord[which(coexp_and_int), "Hub"] <- "Co-expression + Interaction"
    sel_vertex <- c(sel_vertex, coexp_hubs)
  }
  
  colnames(edges) <-  c("X1","Y1","X2","Y2")
  #edges$midX  <- (edges$X1 + edges$X2) / 2
  #edges$midY  <- (edges$Y1 + edges$Y2) / 2
  plotcord[which(plotcord[, "vertex.names"] %in% sel_vertex), "shouldLabel"] <- TRUE
  plotcord$Degree_cut <- cut(plotcord$Degree, breaks=3, labels=FALSE)
  plotcord$in_mod <- TRUE
  #mod_genes <- cem@module[cem@module$modules==name,]$genes
  not_in <- setdiff(plotcord[,"vertex.names"], mod_genes)
  plotcord[which(plotcord[, "vertex.names"] %in% not_in), "in_mod"] <- FALSE
  plotcord$Hub <- ifelse(plotcord$vertex.names %in% MSBB_risk$name,'AD risk gene',plotcord$Hub)
  plotcord$shouldLabel <- ifelse(plotcord$Hub == 'AD risk gene' &
                                   plotcord$in_mod == TRUE,TRUE,plotcord$shouldLabel)
  
  pl <- ggplot(plotcord)  +
    geom_segment(data=edges, aes_(x=~X1, y=~Y1, xend=~X2, yend=~Y2),
                 size = 0.5, alpha=0.5, colour="#DDDDDD") +
    geom_point(aes_(x=~X1, y=~X2, size=~Degree, alpha=~Degree), color=color) +
    geom_label_repel(aes_(x=~X1, y=~X2, label=~vertex.names, color=~Hub),
                     box.padding=unit(1, "lines"),
                     data=function(x){x[x$shouldLabel, ]}) +
    scale_colour_manual(values=c("Co-expression" = "#005E87",
                                 "Interaction" = "#540814",
                                 "Co-expression + Interaction" = "#736E0B",
                                 "AD risk gene" = 'black')) +
    labs(title=name) +
    ggplot2::theme_bw(base_size = 12, base_family = "") +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   legend.key = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(fill = "white",
                                                            colour = NA),
                   panel.border = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
  
  return(pl)
}

modules_colors <- c('#F8766D','#CD9600','#7CAE00',
                    '#00BE67','#00BFC4','#00A9FF',
                    '#C77CFF','#FF61CC')
names(modules_colors) <- mod_names

res <- lapply(mod_names, function(name){
  plot_interaction(ig_obj=human_cem_MSBB@interactions[[name]],
                   n=10, color=modules_colors[name], name=name,
                   mod_genes=module_genes(human_cem_MSBB, name)$genes,
                   coexp_hubs=hubs[[name]],
                   MSBB_risk = AD_risk)
})
names(res) <- mod_names

plots_MSBB <- c('M1','M2','M3','M4')

for(i in 1:length(plots_MSBB)){
  ggsave(file=paste0('NEW_RESULTS/results/cemitool_01/ppi_plots/MSBB/plot_',
                     plots_MSBB[i],'_MSBB_10.jpeg'),
         plot = res[[plots_MSBB[i]]],height = 25,width = 30,units = 'cm',dpi = 600)
}
