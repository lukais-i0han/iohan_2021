#### 
library('dplyr')
library('stringr')
library('RedeR')
library('RColorBrewer')
library('igraph')

mayo_gprofiler_TCX <- read.csv('NEW_RESULTS/gprofiler_results/Human/TCX/mayo_gprofiler_TCX.csv',
                      header = T,
                    stringsAsFactors = F, row.names = 1)
mayo_gprofiler_TCX <- mayo_gprofiler_TCX %>%
  select(query,GO.ID,Genes,Class)



mayo_TCX_edge <- read.csv('NEW_RESULTS/gprofiler_results/mayo_TCX(fdr<0.05)_edge.csv',header = T,
                        stringsAsFactors = F)

mayo_TCX_edge <- mayo_TCX_edge %>% select(EnrichmentMap..Overlap_genes,EnrichmentMap..Overlap_size,
                  name)
colnames(mayo_TCX_edge) <- c('Overlaping_Genes','Overlap_Size','name')

mayo_TCX_edge$name <- str_replace_all(mayo_TCX_edge$name," \\(.*\\) ",'-')

edges_mayo_TCX <- data.frame(str_split_fixed(mayo_TCX_edge$name,'-',2),stringsAsFactors = F)
colnames(edges_mayo_TCX) <- c('node1','node2')

mayo_node_TCX <- read.csv('NEW_RESULTS/gprofiler_results/mayo_TCX(fdr<0.05)_node.csv',header = T, stringsAsFactors = F)
mayo_node_TCX <- mayo_node_TCX %>% select(EnrichmentMap..GS_DESCR,
                            EnrichmentMap..Name)
colnames(mayo_node_TCX) <- c('Description','GO.ID')

mayo_node_TCX <- merge(mayo_node_TCX,mayo_gprofiler_TCX,by='GO.ID')

mayo_node_TCX <- mayo_node_TCX %>% mutate(query = case_when(
  query == 'ALZ[AB.AC.BC]' ~ 'ALZ_ABC',
  query == 'PA[AB.AC.BC]' ~ 'PA_ABC',
  TRUE ~  mayo_node_TCX$query))

mayo_node_TCX$query <- factor(mayo_node_TCX$query,
                  levels = c('ALZ_A','ALZ_B','ALZ_C','ALZ_ABC','PA_A','PA_C',
                      'PA_ABC','PSP_A'))


mayo_node_TCX <- mayo_node_TCX %>% mutate(Group = case_when(
  str_detect(query,'ALZ_') ~ 'ALZ',
  str_detect(query,'PSP_') ~ 'PSP',
  TRUE ~ 'PA'
))


#### Create the igraph object

g <- igraph::graph.data.frame(d=edges_mayo_TCX,directed = F)

### ALZ
rdp <- RedPort() 
calld(rdp)

g_ALZ <- RedeR::subg(g=g, dat = mayo_node_TCX[str_detect(mayo_node_TCX$query,'ALZ_'),],
                      refcol = 1,
                      maincomp = F,connected = T)
g_ALZ <- att.setv(g = g_ALZ, from='Description', to="nodeAlias")
g_ALZ <- att.setv(g = g_ALZ,from="query", to="nodeColor",
                   cols = c('#7FC97F','#FEDF8F','#F0027F','#6990AA'))
g_ALZ <- att.setv(g_ALZ,from='Class',to='nodeShape',shapes = c('ELLIPSE','RECTANGLE','TRIANGLE'))


a<- as_data_frame(g_ALZ,what = 'both')

vertices_ALZ <- a$vertices
edges_ALZ <- a$edges

degree_ALZ <- data.frame(degree(g_ALZ),stringsAsFactors = F)

vertices_ALZ <- vertices_ALZ[rownames(degree_ALZ),]
vertices_ALZ$'degree' <- degree_ALZ$degree.g_ALZ.

h <- igraph::graph.data.frame(d=edges_ALZ,directed = F)



h<-RedeR::subg(g=h, dat = vertices_ALZ,
               refcol = 1,
               maincomp = F,connected = T)

h <- att.setv(h,from='degree',to='nodeSize')

addGraph(rdp,h,theme='tm1')

relax(rdp,p2=400,p5=30,ps=T)

resetd(rdp)


### PSP

g_PSP<- RedeR::subg(g=g, dat = mayo_node_TCX[str_detect(mayo_node_TCX$query,'PSP_'),],
                     refcol = 1,
                     maincomp = F,connected = T)
g_PSP <- att.setv(g = g_PSP, from='Description', to="nodeAlias")
g_PSP <- att.setv(g = g_PSP,from="query", to="nodeColor",
                  cols = c('#7FC97F'))
g_PSP <- att.setv(g_PSP,from='Class',to='nodeShape',shapes = c('ELLIPSE','RECTANGLE','TRIANGLE'))



a<- as_data_frame(g_PSP,what = 'both')

vertices_PSP <- a$vertices
edges_PSP <- a$edges

degree_PSP <- data.frame(degree(g_PSP),stringsAsFactors = F)

vertices_PSP <- vertices_PSP[rownames(degree_PSP),]
vertices_PSP$'degree' <- degree_PSP$degree.g_PSP.

h <- igraph::graph.data.frame(d=edges_PSP,directed = F)

h<-RedeR::subg(g=h, dat = vertices_PSP,
            refcol = 1,
            maincomp = F,connected = T)

h <- att.setv(h,from='degree',to='nodeSize')

addGraph(rdp,h,theme='tm1')

relax(rdp,p2=400,p5=30,ps=T)

resetd(rdp)



relax(rdp,p2=400,p5=30,ps=T)
scl <- g_PSP$legNodeColor$scale 
leg <- g_PSP$legNodeColor$legend 
addLegend.color(rdp, colvec=scl, labvec=leg, 
                title="")

scl <- g_PSP$legNodeShape$shape
leg <- g_PSP$legNodeShape$legend
addLegend.shape(rdp, shapevec = scl, labvec=leg, title="")

resetd(rdp)



#### PA

g_PA<- RedeR::subg(g=g, dat = mayo_node_TCX[str_detect(mayo_node_TCX$query,'PA_'),],
                    refcol = 1,
                    maincomp = F,connected = T)

g_PA <- att.setv(g = g_PA, from='Description', to="nodeAlias")
g_PA <- att.setv(g = g_PA,from="query", to="nodeColor",
                  cols = c('#7FC97F','#FEDF8F','#6990AA'))
g_PA <- att.setv(g_PA,from='Class',to='nodeShape',shapes = c('ELLIPSE','RECTANGLE','TRIANGLE'))


a<- as_data_frame(g_PA,what = 'both')

vertices_PA <- a$vertices
edges_PA <- a$edges

degree_PA <- data.frame(degree(g_PA),stringsAsFactors = F)

vertices_PA <- vertices_PA[rownames(degree_PA),]
vertices_PA$'degree' <- degree_PA$degree.g_PA.

h <- igraph::graph.data.frame(d=edges_PA,directed = F)

h<-RedeR::subg(g=h, dat = vertices_PA,
               refcol = 1,
               maincomp = F,connected = T)

h <- att.setv(h,from='degree',to='nodeSize')

addGraph(rdp,h,theme='tm1')



resetd(rdp)


write.csv(vertices_ALZ,'NEW_RESULTS/results/vertices_ALZ_TCX.csv')
write.csv(vertices_PSP,'NEW_RESULTS/results/vertices_PSP_TCX.csv')
write.csv(vertices_PA,'NEW_RESULTS/results/vertices_PA_TCX.csv')
### RedeR

rdp <- RedPort() 
calld(rdp)

#### CER

mayo_gprofiler_CER <- read.csv('NEW_RESULTS/gprofiler_results/Human/CER/mayo_gprofiler_CER.csv',header = T,
                               stringsAsFactors = F, row.names = 1)
mayo_gprofiler_CER <- mayo_gprofiler_CER %>%
  select(query,GO.ID,Genes,Class)



mayo_CER_edge <- read.csv('NEW_RESULTS/gprofiler_results/mayo_CER(fdr<0.05)_edge.csv',header = T,
                          stringsAsFactors = F)

mayo_CER_edge <- mayo_CER_edge %>% select(EnrichmentMap..Overlap_genes,EnrichmentMap..Overlap_size,
                                          name)
colnames(mayo_CER_edge) <- c('Overlaping_Genes','Overlap_Size','name')

mayo_CER_edge$name <- str_replace_all(mayo_CER_edge$name," \\(.*\\) ",'-')

edges_mayo_CER <- data.frame(str_split_fixed(mayo_CER_edge$name,'-',2),stringsAsFactors = F)
colnames(edges_mayo_CER) <- c('node1','node2')

mayo_node_CER <- read.csv('NEW_RESULTS/gprofiler_results/mayo_CER(fdr<0.05)_node.csv',header = T, stringsAsFactors = F)
mayo_node_CER <- mayo_node_CER %>% select(EnrichmentMap..GS_DESCR,
                                          EnrichmentMap..Name)
colnames(mayo_node_CER) <- c('Description','GO.ID')

mayo_node_CER <- merge(mayo_node_CER,mayo_gprofiler_CER,by='GO.ID')

mayo_node_CER <- mayo_node_CER %>% mutate(query = case_when(
  query == 'ALZ[AB.AC.BC]' ~ 'ALZ_ABC',
  TRUE ~  mayo_node_CER$query))

mayo_node_CER$query <- factor(mayo_node_CER$query,
                              levels = c('ALZ_B','ALZ_C','ALZ_ABC','PA_B','PA_C',
                                         'PSP_A','PSP_B'))
#### Create the igraph object

g <- igraph::graph.data.frame(d=edges_mayo_CER,directed = F)

### ALZ
g_ALZ <- RedeR::subg(g=g, dat = mayo_node_CER[str_detect(mayo_node_CER$query,'ALZ_'),],
                     refcol = 1,
                     maincomp = F,connected = T)
g_ALZ <- att.setv(g = g_ALZ, from='Description', to="nodeAlias")
g_ALZ <- att.setv(g = g_ALZ,from="query", to="nodeColor",
                  cols = c('#FEDF8F','#F0027F','#6990AA'))
g_ALZ <- att.setv(g_ALZ,from='Class',to='nodeShape',shapes = c('ELLIPSE','RECTANGLE','TRIANGLE'))


a<- as_data_frame(g_ALZ,what = 'both')

vertices_ALZ <- a$vertices
edges_ALZ <- a$edges

degree_ALZ <- data.frame(degree(g_ALZ),stringsAsFactors = F)

vertices_ALZ <- vertices_ALZ[rownames(degree_ALZ),]
vertices_ALZ$'degree' <- degree_ALZ$degree.g_ALZ.

h <- igraph::graph.data.frame(d=edges_ALZ,directed = F)

h<-RedeR::subg(g=h, dat = vertices_ALZ,
               refcol = 1,
               maincomp = F,connected = T)

h <- att.setv(h,from='degree',to='nodeSize')

addGraph(rdp,h,theme='tm1')



relax(rdp,p2=400,p5=30,ps=T)
scl <- g_ALZ$legNodeColor$scale 
leg <- g_ALZ$legNodeColor$legend 
addLegend.color(rdp, colvec=scl, labvec=leg, 
                title="")

scl <- g_ALZ$legNodeShape$shape
leg <- g_ALZ$legNodeShape$legend
addLegend.shape(rdp, shapevec = scl, labvec=leg, title="")

resetd(rdp)


### PSP

g_PSP<- RedeR::subg(g=g, dat = mayo_node_CER[str_detect(mayo_node_CER$query,'PSP_'),],
                    refcol = 1,
                    maincomp = F,connected = T)
g_PSP <- att.setv(g = g_PSP, from='Description', to="nodeAlias")
g_PSP <- att.setv(g = g_PSP,from="query", to="nodeColor",
                  cols = c('#7FC97F','#F0027F'))
g_PSP <- att.setv(g_PSP,from='Class',to='nodeShape',shapes = c('ELLIPSE','RECTANGLE','TRIANGLE'))



a<- as_data_frame(g_PSP,what = 'both')

vertices_PSP <- a$vertices
edges_PSP <- a$edges

degree_PSP <- data.frame(degree(g_PSP),stringsAsFactors = F)

vertices_PSP <- vertices_PSP[rownames(degree_PSP),]
vertices_PSP$'degree' <- degree_PSP$degree.g_PSP.

h <- igraph::graph.data.frame(d=edges_PSP,directed = F)

h<-RedeR::subg(g=h, dat = vertices_PSP,
               refcol = 1,
               maincomp = F,connected = T)

h <- att.setv(h,from='degree',to='nodeSize')

addGraph(rdp,h,theme='tm1')


relax(rdp,p2=400,p5=30,ps=T)
scl <- g_PSP$legNodeColor$scale 
leg <- g_PSP$legNodeColor$legend 
addLegend.color(rdp, colvec=scl, labvec=leg, 
                title="")

scl <- g_PSP$legNodeShape$shape
leg <- g_PSP$legNodeShape$legend
addLegend.shape(rdp, shapevec = scl, labvec=leg, title="")

resetd(rdp)



#### PA

g_PA<- RedeR::subg(g=g, dat = mayo_node_CER[str_detect(mayo_node_CER$query,'PA_'),],
                   refcol = 1,
                   maincomp = F,connected = T)
g_PA <- att.setv(g = g_PA, from='Description', to="nodeAlias")
g_PA <- att.setv(g = g_PA,from="query", to="nodeColor",
                 cols = c('#F0027F','#6990AA'))
g_PA <- att.setv(g_PA,from='Class',to='nodeShape',shapes = c('ELLIPSE','RECTANGLE','TRIANGLE'))

a<- as_data_frame(g_PA,what = 'both')

vertices_PA <- a$vertices
edges_PA <- a$edges

degree_PA <- data.frame(degree(g_PA),stringsAsFactors = F)

vertices_PA <- vertices_PA[rownames(degree_PA),]
vertices_PA$'degree' <- degree_PA$degree.g_PA.

h <- igraph::graph.data.frame(d=edges_PA,directed = F)

h<-RedeR::subg(g=h, dat = vertices_PA,
               refcol = 1,
               maincomp = F,connected = T)

h <- att.setv(h,from='degree',to='nodeSize')

addGraph(rdp,h,theme='tm1')


relax(rdp,p2=400,p5=30,ps=T)
scl <- g_PA$legNodeColor$scale 
leg <- g_PA$legNodeColor$legend 
addLegend.color(rdp, colvec=scl, labvec=leg, 
                title="")

scl <- g_PA$legNodeShape$shape
leg <- g_PA$legNodeShape$legend
addLegend.shape(rdp, shapevec = scl, labvec=leg, title="")

resetd(rdp)

write.csv(vertices_ALZ,'NEW_RESULTS/results/vertice_ALZ_CER.csv')
write.csv(vertices_PSP,'NEW_RESULTS/results/vertice_PSP_CER.csv')
write.csv(vertices_PA,'NEW_RESULTS/results/vertice_PA_CER.csv')
