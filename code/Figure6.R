rm(list = ls())
library(igraph)
library(qgraph)
library(reshape2)
library(grid)

load("SAD8_module.Rdata") #result from ODE_solving

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
  
  net2 <- graph_from_data_frame( d=after,vertices = nodes,directed = T )
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
               vertex.label.cex=V(net)$ind_effect*2.2,
               vertex.size=V(net)$ind_effect*30+5,
               edge.arrow.size=E(net)$Effect*3,
               edge.curved=0.05,
               edge.color=E(net)$edge.colour,
               edge.frame.color=E(net)$edge.colour,
               edge.width=E(net)$Effect*10,
               vertex.color=V(net)$node.colour,
               layout=l,
               main=title,
               margin=c(-.05,-.05,-.05,-.05)
  )
  box("figure")
  return(net2)
}


get_summary <- function(net_all){
  
  
  after_n = data.frame(do.call(rbind, lapply(net_all, "[[", "edge")))
  
  in_n = data.frame(table(after_n[,2]))
  out_n = data.frame(table(after_n[,1]))
  
  #levels(in_n$Var1) = in_n$Var1[order(in_n$Freq)]
  #levels(in_n$Var1) = order(in_n$Freq)
  df = Reduce(function(x, y) merge(x, y, all=TRUE, by ="Var1"), list(in_n, out_n))
  colnames(df) = c("cluster","in1","out1")
  df[is.na(df)] <- 0
  df[,3] <- -df[,3]
  
  
  g1 <- ggplot() +
    geom_col(df, mapping=aes(reorder(cluster, out1), out1),
             position = position_dodge(), fill = "#65C18C",width = 1, size = 0.3, alpha = 0.67, show.legend = T) +
    geom_col(df, mapping=aes(cluster, in1),
             position = position_dodge(), fill = "#FF7BA9",width = 1, size = 0.3, alpha = 0.67, show.legend = T) +
    theme(panel.grid = element_blank(), 
          panel.background = element_rect(color = 'black', fill = 'transparent')) +
    xlab(NULL) + ylab(NULL)+ 
    scale_color_manual(name=NULL,
                       breaks=c('Incoming', 'Outgoing'),
                       values = c("Incoming" = "#65C18C", "Outgoing" = "#FF7BA9")) +
    #annotate('text', label = 'Outgoing', 9, 13, size = 6) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12)) +
    theme(legend.position = c(0.9, 0.9)) + 
    scale_y_continuous(limits = c(-35,35), breaks = c(-30,-20,-10, 0, 10, 20,30),labels = c(30,20,10,0,10,20,30)) 
  #geom_hline(yintercept = 0, size = 0.5)
  
  return(g1)
}

get_summary2 <- function(net, j ){
  a = eigen_centrality(net,
                       directed = TRUE,
                       weights = abs(E(net)$Effect),
                       scale = FALSE)
  name = c("STR","TSS","TTR","TTS","SSR","SSS","STS","TSR")
  d = data.frame(names = names(a$vector), values = a$vector)
  d = d[-(which(d$values < 0.01)),]
  p1 = ggplot(d) + geom_col(mapping=aes(reorder(names,  values, decreasing = T),  values))+ theme_bw()+
    ylab("Centrity") + xlab("Nodes") +
    scale_y_continuous(expand = c(0,NA))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12))+
    ggtitle(name[j]) + theme(plot.title = element_text(hjust = 0.5))
  p1
}

get_summary3 <- function(net, j ){
  net_all = lapply(1:8,function(c)network_plot(Module_net[[c]]))
  
  a = lapply(c(6,5,4,3,7,8,2,1),function(c) eigen_centrality(net_all[[c]],
                                                             directed = TRUE,
                                                             weights = abs(E(net_all[[c]])$Effect),
                                                             scale = FALSE))
  name = c("SSS","SSR","TTS","TTR","STS","TSR","TSS","STR")
  names(a) = name
  
  
  d = lapply(1:8,function(x)data.frame(names = names(a[[x]]$vector), 
                                       values = a[[x]]$vector))
  d = Reduce(cbind,d)
  d2 = d[,seq(2,16,2)]
  colnames(d2) = name
  
  d2$STS = -d2$STS
  d2$TSR = -d2$TSR
  d2$TSS = -d2$TSS
  d2$STR = -d2$STR
  
  
  d3 = d2[which(apply(d2>0.05,1,any)==T),]
  #d3$ID = rownames(d3)
  d4 = melt(as.matrix(d3))
  d4$group = c(rep("up",28),rep("down",28))
  d4$group2 = c(rep(1,7),rep(2,7),rep(3,7),rep(4,7),
                rep(1,7),rep(2,7),rep(3,7),rep(4,7))
  
  
  d5 = list(d4[c(1:7,29:35),],
            d4[c(8:14,36:42),],
            d4[c(15:21,43:49),],
            d4[c(22:28,50:56),])
  
  d = d5[[1]]
  d$type = substr(as.character(d$Var2),1,2)
  
  all_name1 = c("S/S","S/S","T/T","T/T")
  all_name2 = c("S/T","T/S","T/S","S/T")
  
  cols <- c("ST" = '#FFB534', "TS" = '#65B741', "SS" = "#FF004D", "TT" = "#3652AD")
  
  plot_p <- function(i){
    d = d5[[i]]
    d$type = substr(as.character(d$Var2),1,2)
    p = ggplot(data = d, mapping = aes(reorder(Var1,  abs(value), decreasing = T),value)) + 
      geom_bar(mapping = aes(fill = type),stat = 'identity')+
      scale_fill_manual(values = cols)+
      theme_bw()+ xlab(NULL)+ylab(NULL)+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12))+
      scale_y_continuous(limits = c(-1,1),breaks = c(-1,-0.5,0,0.5,1))+ 
      theme(legend.position="none")+
      annotate("text", x = 3.5, y = 0.9, label = all_name1[i])+
      annotate("text", x = 3.5, y = -0.9, label = all_name2[i])
    p
  }
  
  #plot_p(1)
  p = lapply(1:4,function(x)plot_p(x))
  
  pp = wrap_plots(p, ncol = 4)
  
  pp
}



ggsave("Fig6B.pdf",pp,width = 18, height = 4)



pdf("Fig6A.pdf",width=20,height=15)
#vp1 <- viewport(height=unit(0.3, "npc"), width=unit(0.35, "npc"), y=0.15, x=0.25)
#vp2 <- viewport(height=unit(0.3, "npc"), width=unit(0.35, "npc"), y=0.15, x=0.75)
layout(matrix(1:12, 3, 4, byrow = TRUE))
network_plot(Module_net[[6]], title = "SSS")
network_plot(Module_net[[4]], title = "TTS")
network_plot(Module_net[[7]], title = "STS")
network_plot(Module_net[[2]], title = "TSS")

network_plot(Module_net[[5]], title = "SSR")
network_plot(Module_net[[3]], title = "TTR")
network_plot(Module_net[[1]], title = "STR")
network_plot(Module_net[[8]], title = "TSR")
#print(p1, vp=vp1)
#print(p2, vp=vp2)
dev.off()
