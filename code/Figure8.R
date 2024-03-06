rm(list = ls())
library(topGO)
library(Biostrings)
library(pbapply)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(grid)
library(reshape2)
library(RColorBrewer)
library(rCharts)
#from previous calculation
load(file = "SAD8_subcluster1.Rdata")
load(file = "SAD8_subcluster2.Rdata")

#result = result_31_40[[9]]

GOenrich <-  function(deglist,geneID2GO,useFDR = FALSE,cut=0.05, outdir = NULL,outprefix=NUL) {
  if (is.null(outdir)) {
    outdir = getwd()
  }
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
  BP.res = GOenrichsub(deglist,
                       geneID2GO,
                       'BP',
                       paste0(outdir, '/', outprefix))
  MF.res = GOenrichsub(deglist,
                       geneID2GO,
                       'MF',
                       paste0(outdir, '/', outprefix))
  CC.res = GOenrichsub(deglist,
                       geneID2GO,
                       'CC',
                       paste0(outdir, '/', outprefix))
  
  all = rbind(BP.res, MF.res, CC.res)
  all[all$Pvalue=='< 1e-30','Pvalue']=1e-31
  sigres=all[all$Pvalue < cut,]
  return(list(all,sigres))
}

GOenrichsub <- function(deglist, geneID2GO, class, output) {
  cat('Processing', class, 'class...')
  geneList <- rep(1, length(geneID2GO))
  names(geneList) <- names(geneID2GO)
  geneList[match(deglist, names(geneList))] = 0
  sampleGOdata <-
    suppressMessages(
      new(
        "topGOdata",
        nodeSize = 1,
        ontology = class,
        allGenes = geneList,
        annot = annFUN.gene2GO,
        gene2GO = geneID2GO,
        geneSel = function(allScore) {
          return(allScore < 0.01)
        }
      )
    )
  cat('...')
  result <-
    suppressMessages(runTest(
      sampleGOdata,
      algorithm = "elim",
      statistic = "fisher"
    ))
  allRes <-
    GenTable(
      sampleGOdata,
      Pvalue = result,
      orderBy = "Pvalue",
      topNodes = attributes(result)$geneData[4],
      numChar = 1000
    )
  cat('...')
  
  gene = NULL
  GOs = sampleGOdata@graph@nodeData@data
  for (i in allRes$GO.ID) {
    e = paste0('g=names(GOs$`', i, '`$genes)')
    eval(parse(text = e))
    
    d = intersect(g, deglist)
    e = paste0('gene[\'', i, '\']=toString(d)')
    eval(parse(text = e))
  }
  allRes$gene = gene
  allRes = cbind(allRes,
                 AllDEG = length(geneList[geneList == 0]),
                 Class = class)
  cat('done!\n')
  return(allRes)
}

get_db <- function(){
  #dat = read.delim("GO.out", header = T, sep="\t")
  #GOcontent = data.frame(cbind(GID=dat$qpid,GOID=paste0('GO:', sprintf("%07d", dat$goid))))
  #geneID2GO =lapply(unique(GOcontent$GID), function (x) GOcontent$GOID[GOcontent$GID==x])
  #names(geneID2GO) = unique(GOcontent$GID)
  dat = read.delim("contig_p94_MBY_GO.txt", header = T, sep = "\t")
  get_list <- function(i){
    tmp = strsplit(dat[i,]$GO.IDs,"; ")
    c(unlist(lapply(tmp, substr,3,12)))
  }
  GOcontent = lapply(1:nrow(dat),function(c) get_list(c))
  #GOcontent = lapply(1:10,function(c) get_list(c))
  names(GOcontent) = unique(dat$SeqName)
  return(GOcontent)
}

geneID2GO = get_db()

#dat = read.delim(file = "contig_p94_MBY_GO.txt",sep = "\t")

get_gene_list <- function(result,i){
  df_cluster = result$cluster[which(result$cluster$`apply(omega, 1, which.max)`==i),]
  df_cluster = df_cluster[,-65]
  return(rownames(df_cluster))
}



genelist = list(get_gene_list(result_M3[[9]], 6),
                get_gene_list(result_M5[[11]], 2),
                get_gene_list(result_M7[[16]], 1),
                get_gene_list(result_M8[[10]], 2),
                get_gene_list(result_M13[[12]], 9),
                get_gene_list(result_M15[[9]], 4),
                get_gene_list(result_M16[[19]], 17),
                get_gene_list(result_M17[[19]], 4),
                get_gene_list(result_M18[[9]], 6),
                get_gene_list(result_M19[[7]], 6),
                get_gene_list(result_M21[[16]], 2),
                get_gene_list(result_M22[[13]], 2),
                get_gene_list(result_M24[[17]], 8),
                get_gene_list(result_M25[[7]],51),
                get_gene_list(result_M27[[16]], 1),
                get_gene_list(result_M28[[8]], 3),
                get_gene_list(result_M29[[17]], 2),
                get_gene_list(result_M30[[12]], 1),
                get_gene_list(result_M32[[11]], 1),
                get_gene_list(result_M36[[19]], 16) )
#genelist = dat$SeqName[1:1000]
#save.image(file = "GO.RData")


#GO_a <- GOenrich(genelist[[7]],geneID2GO)

GO_all <- lapply(1:length(genelist),function(c) GOenrich(genelist[[c]],geneID2GO))


#dataset = GO_a
#####plot a module
GObubble <- function(dataset,cut = 0.05,top = 10, MainTitle = NULL) {
  require(tidyverse)
  
  dataset$Pvalue = -log10(as.numeric(dataset$Pvalue))
  dataset = dataset %>% group_by(Class) %>% slice_max(order_by = Pvalue, n = top) %>% as.data.frame
  
  dataset$GeneRatio = dataset$Significant / dataset$AllDEG * 100
  order <- order(dataset$GeneRatio)
  dataset$Term <- factor(dataset$Term, levels = dataset$Term[order])
  xLabel = 'Count'
  yLabel = ''
  sizeLabel = 'Gene Ratio (%)'
  colorLabel = expression('-log'[10] * 'P-value')
  DotRange = c(2, 10)
  col = colorRampPalette(c("#B4FF9F","#F9F3DF", "#FFA1A1"))(100)
  col = col[15:85]
  
  colorbreaks = pretty(dataset$Pvalue, 5)
  sizebreaks = pretty(dataset$GeneRatio, 5)
  
  g <- ggplot(data = dataset, mapping = aes(Significant, Term)) +
    geom_point(aes(size = GeneRatio, color = Pvalue),
               show.legend = T,
               stroke = 1) +
    facet_grid(Class ~ ., scales = "free", space = "free") +
    scale_size(
      name = sizeLabel,
      breaks = sizebreaks,
      labels = sizebreaks,
      range = DotRange,
      guide = "legend"
    ) +
    scale_color_gradientn(colours = col, breaks = colorbreaks) +
    scale_y_discrete(position = "left") +
    labs(
      title = MainTitle,
      x = xLabel,
      y = yLabel,
      size = sizeLabel,
      color = colorLabel
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(colour = "black", size = 11),
      axis.text.y = element_text(lineheight = 0.65),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.position = "right",
      strip.background = element_blank(),
      strip.text = element_text(size = 15, face = 'bold')
    )
  
  return(g)
}


name_all = c("SM3_6","SM5_2","SM7_1","SM8_2","SM13_9","SM15_4","SM16_17",
             "SM17_4","SM18_9","SM19_6","SM21_2","SM22_2","SM24_8",
             "SM25_5","SM27_1","SM28_3","SM29_2","SM30_1","SM32_1","SM36_16")

p_list = lapply(1:length(GO_all), function(c) GObubble(GO_all[[c]][[2]], MainTitle = name_all[c]))

pp = wrap_plots(p_list[1:2],nrow = 2)
ggsave("g1.pdf",pp ,width = 16, height = 13)

pp = wrap_plots(p_list[13:20],ncol = 2)
ggsave("g3.pdf",pp ,width = 22, height = 20)


##########
                normalization <- function(x, z = 0.2){(x-min(x))/(max(x)-min(x))+z}

darken <- function(color, factor=1.2){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

network_conversion <- function(result){
  n = ncol(result$fit)
  
  effect.mean = apply(result$fit,2,mean)[4:n]
  effect.predict.mean = apply(result$predict,2,mean)[4:n]
  effect.total = colSums(result$fit)[4:n]
  
  temp = matrix(NA,nrow = n-3, ncol=3)
  colnames(temp) = c("From", "To", "Effect")
  temp[,1] = colnames(result$fit)[4:n]
  temp[,2] = colnames(result$fit)[4]
  temp[,3] = effect.predict.mean
  if (nrow(temp)==2) {
    temp = t(data.frame(temp[-1,]))
  } else{
    temp = data.frame(temp[-1,])
  }
  
  temp2 = matrix(NA,nrow = n-3, ncol=3)
  colnames(temp2) = c("From", "To", "Effect")
  temp2[,1] = colnames(result$fit)[4:n]
  temp2[,2] = colnames(result$fit)[4]
  temp2[,3] = effect.total
  if (nrow(temp2)==2) {
    temp2 = t(data.frame(temp2[-1,]))
  } else{
    temp2 = data.frame(temp2[-1,])
  }
  output <- list(ind.name = colnames(result$fit)[4],
                 dep.name = colnames(result$fit)[5:n],
                 ODE.par = result$ODE.value,
                 ind.par = result$LOP_par[,3],
                 dep.par = result$LOP_par[,4:(n-1)],
                 effect.mean = effect.predict.mean,
                 effect.total = effect.total,
                 effect.all = result$fit,
                 edge = temp,
                 edge.total = temp2,
                 ind.effect = effect.predict.mean[1])
  return(output)
}

#result = Module_net[[1]]
network_plot <- function(result, title = NULL, maxeffect = NULL, type = NULL){
  source("igraphplot2.R")
  environment(plot.igraph2) <- asNamespace('igraph')
  environment(igraph.Arrows2) <- asNamespace('igraph')
  extra <- sapply(result,"[[", "ind.effect")
  
  if (is.null(type)) {
    after <- data.frame(do.call(rbind, lapply(result, "[[", "edge")))
  } else{
    after <- data.frame(do.call(rbind, lapply(result, "[[", "edge.total")))
    
  }
  after$Effect = as.numeric(after$Effect)
  rownames(after) = NULL
  
  
  after$edge.colour = NA
  for (i in 1:nrow(after)) {
    if(after$Effect[i]>=0){
      after$edge.colour[i] = "#FE433C"
    } else{
      after$edge.colour[i] = "#0095EF"
    }
  }
  
  #nodes
  nodes <- data.frame(unique(after[,2]),unique(after[,2]),extra)
  colnames(nodes) <- c("id","name","ind_effect")
  nodes$influence <- aggregate(Effect ~ To, data = after, sum)[,2]
  nodes$node.colour = NA
  for (i in 1:nrow(nodes)) {
    if(nodes$influence[i]>=0){
      nodes$node.colour[i] = "#FFC4C4"
    } else{
      nodes$node.colour[i] = "#89CFFD"
    }
  }
  
  #normalization
  if (is.null(maxeffect)) {
    after[,3] <- normalization(abs(after[,3]))
    nodes[,3:4] <- normalization(abs(nodes[,3:4]))
  } else{
    after[,3] <- (abs(after[,3]))/maxeffect*1.5+0.1
    nodes[,3:4] <- (abs(nodes[,3:4]))/maxeffect*1.5+0.1
  }
  
  net <- graph_from_data_frame( d=after,vertices = nodes,directed = T )
  
  #layout
  set.seed(2)
  g <- barabasi.game(nrow(nodes), directed=FALSE)
  e <- get.edgelist(g,names=FALSE)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),
                                         area=8*(vcount(g)^2),repulse.rad=(vcount(g)^3.1))
  
  #l <- layout_randomly(net)
  #l = layout_as_star(net)
  plot.igraph2(net,
               vertex.label=V(net)$name,
               vertex.label.color="black",
               vertex.shape="circle",
               vertex.frame.color=NA,
               vertex.label.cex=V(net)$ind_effect*2+0.2,
               vertex.size=V(net)$ind_effect*20+5,
               edge.arrow.size=E(net)$Effect*1,
               edge.curved=0.05,
               edge.color=E(net)$edge.colour,
               edge.frame.color=E(net)$edge.colour,
               edge.width=E(net)$Effect*3,
               vertex.color=V(net)$node.colour,
               layout=l,
               main=title,
               margin=c(-.05,-.05,-.05,-.05)
  )
  net
  #box("figure")
}
Module_net = lapply(1:8,function(c) lapply(gene_ode_submodule3[[6]][[c]]$ode_result, network_conversion))
net = network_plot(Module_net[[1]])
                
edgelist <- get.data.frame(net) 
colnames(edgelist) <- c("source","target","value")

edgelist$source <- as.character(edgelist$source)
edgelist$target <- as.character(edgelist$target)

sankeyPlot <- rCharts$new()
sankeyPlot$setLib('libraries/widgets/d3_sankey')
sankeyPlot$setTemplate(script = "libraries/widgets/d3_sankey/layouts/chart.html")

sankeyPlot$set(
  data = edgelist,
  nodeWidth = 15,
  nodePadding = 10,
  layout = 32,
  width = 960,
  height = 500
)
