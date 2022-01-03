library('dplyr')
library('stringr')
library('RedeR')
library('RColorBrewer')
library('igraph')

FAD <- read.csv('NEW_RESULTS/gprofiler_results/5XFAD/fad5x_gprofiler2.csv',
                               header = T,
                               stringsAsFactors = F,row.names = 1)
FAD <- FAD %>%
  select(query,GO.ID,Genes,Class)



FAD_edge <- read.csv('NEW_RESULTS/gprofiler_results/5XFAD/FAD_enrichment_edge.csv',header = T,
                          stringsAsFactors = F)

FAD_edge <- FAD_edge %>% select(EnrichmentMap..Overlap_genes,EnrichmentMap..Overlap_size,
                                          name)

colnames(FAD_edge) <- c('Overlaping_Genes','Overlap_Size','name')

FAD_edge$name <- str_replace_all(FAD_edge$name," \\(.*\\) ",'-')


FAD_edge <- data.frame(str_split_fixed(FAD_edge$name,'-',2),stringsAsFactors = F)
colnames(FAD_edge) <- c('node1','node2')

FAD_node <- read.csv('NEW_RESULTS/gprofiler_results/5XFAD/CCX_enrichment_node.csv',header = T, stringsAsFactors = F)
FAD_node <- FAD_node %>% select(EnrichmentMap..GS_DESCR,
                                          EnrichmentMap..Name)
colnames(FAD_node) <- c('Description','GO.ID')

FAD_node <- merge(FAD_node,FAD,by='GO.ID')


FAD_node <- FAD_node[FAD_node$GO.ID %in% mayo_node_TCX$GO.ID,]
#### Create the igraph object

g <- igraph::graph.data.frame(d=FAD_edge,directed = F)

rdp <- RedPort() 
calld(rdp)



### CCX
g_CCX <- RedeR::subg(g=g, dat = FAD_node[str_detect(FAD_node$query,'CCX'),],
                     refcol = 1,
                     maincomp = T,connected = T)
g_CCX <- att.setv(g = g_CCX, from='Description', to="nodeAlias")
g_CCX <- att.setv(g = g_CCX,from="query", to="nodeColor",
                  cols = c('#F0027F','#6990AA','#FEDF8F','#7FC97F'))
g_CCX <- att.setv(g_CCX,from='Class',to='nodeShape',
                shapes = c('ELLIPSE','RECTANGLE','TRIANGLE'))


a<- as_data_frame(g_CCX,what = 'both')

vertices_CCX <- a$vertices
edges_CCX <- a$edges

degree_CCX <- data.frame(degree(g_CCX),stringsAsFactors = F)

vertices_CCX <- vertices_CCX[rownames(degree_CCX),]
vertices_CCX$'degree' <- degree_CCX$degree.g_CCX.

h <- igraph::graph.data.frame(d=edges_CCX,directed = F)


h<-RedeR::subg(g=h, dat = vertices_CCX,
               refcol = 1,
               maincomp = F,connected = T)

h <- att.setv(h,from='degree',to='nodeSize')

addGraph(rdp,h,theme='tm1')


relax(rdp,p2=400,p5=30,ps=T)
scl <- g_CCX$legNodeColor$scale 
leg <- g_CCX$legNodeColor$legend 
addLegend.color(rdp, colvec=scl, labvec=leg, 
                title="")

scl <- g_CCX$legNodeShape$shape
leg <- c(g_CCX$legNodeShape$legend)
addLegend.shape(rdp, shapevec = scl, labvec=leg, title="")

resetd(rdp)

### HIP

g_HIP <- RedeR::subg(g=g, dat = FAD_node[str_detect(FAD_node$query,'HIP'),],
                     refcol = 1,
                     maincomp = F,connected = F)
g_HIP <- att.setv(g = g_HIP, from='Description', to="nodeAlias")
g_HIP <- att.setv(g = g_HIP,from="query", to="nodeColor",
                  cols = c('#F0027F','#6990AA','#FEDF8F','#7FC97F'))
g_HIP <- att.setv(g_HIP,from='Class',to='nodeShape',shapes = c('ELLIPSE','RECTANGLE','TRIANGLE'))



a<- as_data_frame(g_HIP,what = 'both')

vertices_HIP <- a$vertices
edges_HIP <- a$edges

degree_HIP <- data.frame(degree(g_HIP),stringsAsFactors = F)

vertices_HIP <- vertices_HIP[rownames(degree_HIP),]
vertices_HIP$'degree' <- degree_HIP$degree.g_HIP.

h <- igraph::graph.data.frame(d=edges_HIP,directed = F)

h<-RedeR::subg(g=h, dat = vertices_HIP,
               refcol = 1,
               maincomp = F,connected = T)

h <- att.setv(h,from='degree',to='nodeSize')

addGraph(rdp,h,theme='tm1')


relax(rdp,p2=400,p5=30,ps=T)
scl <- g_HIP$legNodeColor$scale 
leg <- g_HIP$legNodeColor$legend 
addLegend.color(rdp, colvec=scl, labvec=leg, 
                title="")

scl <- g_HIP$legNodeShape$shape
leg <- c(g_HIP$legNodeShape$legend)
addLegend.shape(rdp, shapevec = scl, labvec=leg, title="")

resetd(rdp)



#### Taud35
TAU <- read.csv('NEW_RESULTS/gprofiler_results/TAUd35/taud35_gprofiler.csv',
                header = T,
                stringsAsFactors = F,row.names = 1)

TAU <- TAU %>%
  select(query,GO.ID,Genes,Class)



TAU_edge <- read.csv('NEW_RESULTS/gprofiler_results/TAUd35/TAU_enrichment_edge.csv',header = T,
                     stringsAsFactors = F)

TAU_edge <- TAU_edge %>% select(EnrichmentMap..Overlap_genes,EnrichmentMap..Overlap_size,
                                name)
colnames(TAU_edge) <- c('Overlaping_Genes','Overlap_Size','name')

TAU_edge$name <- str_replace_all(TAU_edge$name," \\(.*\\) ",'-')

TAU_edge <- data.frame(str_split_fixed(TAU_edge$name,'-',2),stringsAsFactors = F)
colnames(TAU_edge) <- c('node1','node2')

TAU_node <- read.csv('NEW_RESULTS/gprofiler_results/TAUd35/TAU_enrchiment_node.csv',header = T, stringsAsFactors = F)
TAU_node <- TAU_node %>% select(EnrichmentMap..GS_DESCR,
                                EnrichmentMap..Name)
colnames(TAU_node) <- c('Description','GO.ID')

TAU_node <- merge(TAU_node,TAU,by='GO.ID')
TAU_node <- TAU_node[str_detect(TAU_node$Description,'syna|neur|channel|immune'),]

#### Create the igraph object

g <- igraph::graph.data.frame(d=TAU_edge,directed = F)

g_TAU <- RedeR::subg(g=g, dat = TAU_node,
                     refcol = 1,
                     maincomp = F,connected = F)
g_TAU <- att.setv(g = g_TAU, from='Description', to="nodeAlias")
g_TAU <- att.setv(g = g_TAU,from="query", to="nodeColor",
                  cols = c('#6990AA','#7FC97F'))
g_TAU <- att.setv(g_TAU,from='Class',to='nodeShape',shapes = c('ELLIPSE','RECTANGLE','TRIANGLE'))



a<- as_data_frame(g_TAU,what = 'both')

vertices_TAU <- a$vertices
edges_TAU <- a$edges

degree_TAU <- data.frame(degree(g_TAU),stringsAsFactors = F)

vertices_TAU <- vertices_TAU[rownames(degree_TAU),]
vertices_TAU$'degree' <- degree_TAU$degree.g_TAU.

h <- igraph::graph.data.frame(d=edges_TAU,directed = F)

h<-RedeR::subg(g=h, dat = vertices_TAU,
               refcol = 1,
               maincomp = F,connected = T)

h <- att.setv(h,from='degree',to='nodeSize')


addGraph(rdp,h,theme='tm1')


relax(rdp,p2=400,p5=30,ps=T)
scl <- g_TAU$legNodeColor$scale 
leg <- g_TAU$legNodeColor$legend 
addLegend.color(rdp, colvec=scl, labvec=leg, 
                title="")

scl <- g_TAU$legNodeShape$shape
leg <- c(g_TAU$legNodeShape$legend)
addLegend.shape(rdp, shapevec = scl, labvec=leg, title="")

resetd(rdp)


write.csv(vertices_CCX,'NEW_RESULTS/results/vertices_CCX.csv')
write.csv(vertices_HIP,'NEW_RESULTS/results/vertices_HIP.csv')
write.csv(vertices_TAU,'NEW_RESULTS/results/vertices_TAU.csv')
