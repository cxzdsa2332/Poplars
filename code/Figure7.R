rm(list = ls())
library(igraph)
library(qgraph)
library(reshape2)
library(grid)
library(ggplot2)
library(ggrepel)
library(patchwork)

#load("network/sad8/SAD8_module.Rdata") #result from ODE solving

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


Module_net = lapply(1:8,function(c) lapply(gene_ode_submodule3[[6]][[c]]$ode_result, network_conversion))

p1 = get_summary(Module_net[[1]])


n1 = network_plot(Module_net[[1]], title = "STR")
get_summary2 <- function(net){
  a = eigen_centrality(net,directed = FALSE,scale = TRUE,
                       #weights = V(net)$ind_effect,
                       options = arpack_defaults)
  d = data.frame(names = names(a$vector), values = a$vector)
  #p1 = ggplot(d) + geom_col(mapping=aes(reorder(names,  values, decreasing = T),  values))+ theme_bw()+
  #  xlab("Centrity") + ylab("Nodes") +
  #  scale_y_continuous(expand = c(0,NA))+
  #  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12))
  #p1
  d
}

d = get_summary2(n1)

key_gene = d[which(d$values>0.9),]
result = Module_net[[1]]
network_plot1 <- function(result, title = NULL, maxeffect = NULL, type = NULL){
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
    after[,3] <- (abs(after[,3]))/maxeffect*1.5+0.05
    nodes[,3:4] <- (abs(nodes[,3:4]))/maxeffect*1.5+0.1
  }
  
  
  
  #layout
  set.seed(2)
  g <- barabasi.game(nrow(nodes), directed=FALSE)
  e <- get.edgelist(g,names=FALSE)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),
                                         area=8*(vcount(g)^2),repulse.rad=(vcount(g)^3.1))
  
  #added
  nodes$name = NA
  nodes$name[match(key_gene$names,nodes$id)] = nodes$id[match(key_gene$names,nodes$id)]
  nodes$node.colour[match(key_gene$names,nodes$id)] = "#FFED00"
  nodes$board = NA
  nodes$board[match(key_gene$names,nodes$id)] = "black"
  #l <- layout_randomly(net)
  #l = layout_as_star(net)
  net <- graph_from_data_frame( d=after,vertices = nodes,directed = T )
  plot.igraph2(net,
               vertex.label=V(net)$name,
               #vertex.label=NA,
               vertex.frame.color=nodes$board,
               vertex.label.color="black",
               vertex.shape="circle",
               vertex.label.cex=V(net)$ind_effect*3,
               vertex.size=V(net)$ind_effect*8+1,
               edge.arrow.size=E(net)$Effect*0.15,
               edge.curved=0.05,
               edge.color=E(net)$edge.colour,
               edge.frame.color=E(net)$edge.colour,
               edge.width=E(net)$Effect,
               vertex.color=V(net)$node.colour,
               layout=l,
               main=title,
               margin=c(-.05,-.05,-.05,-.05)
  )
  #box("figure")
}

#network_plot1(Module_net[[1]], title = "STR")

#result = Module_net

Module_net = lapply(1:8,function(c) gene_ode_submodule3[[6]][[c]]$ode_result)
select = 2
gene = lapply(1:8,function(c)
  Module_net[[c]][[which(names(Module_net[[c]])==key_gene$names[select])]])

result = gene[[1]]
qdODEplot_convert <- function(result){
  data = result$predict
  n = ncol(data)
  colnames(data)[4:n] = c(paste0("ind.",colnames(data)[4]),
                          paste0("dep.",colnames(data)[5:n]))
  
  plot.df = melt(data, id.vars = c("x"))
  
  name = levels(plot.df[,2])
  
  ind.name = name[grep("ind", name)]
  ind.name2 = strsplit(ind.name,split = "\\.")[[1]][2]
  ind.df <- subset(plot.df, plot.df[,2] == ind.name)
  ind.df$type = "ind"
  ind.df$variable = ind.name2
  
  depname = levels(plot.df[,2])[grep("dep",name )]
  dep.df <- subset(plot.df, plot.df[,2] %in% depname)
  dep.df$type = "dep"
  dep.df$variable = sapply(strsplit(as.character(dep.df$variable),"\\."),"[",2)
  
  
  original.df = subset(plot.df, plot.df[,2] == "y")
  original.df$type = "original"
  
  fit.df = subset(plot.df, plot.df[,2] == "y.fit")
  fit.df$type = "fit"
  
  plot.df2 = rbind(ind.df, dep.df,fit.df)
  
  name.df = subset(plot.df2, plot.df[,1] == max(plot.df2[,1]))
  name.df = name.df[-nrow(name.df),]
  name.df[,2][name.df[,2] == "y.fit"] = ind.name2
  
  name.df = name.df[-which(name.df[,4] == "fit"),]
  
  name.df[,1] = name.df[,1]*1.00
  return_obj = list(plot.df2 = plot.df2,
                    name.df = name.df,
                    ind.name2 = ind.name2)
  return(return_obj)
}

qdODE_plot_base <- function(result,label = 10, show.legend = TRUE){
  result2 = qdODEplot_convert(result)
  plot.df2 = result2$plot.df2
  name.df = result2$name.df
  ind.name2 = result2$ind.name2
  name.df = name.df[-1,]
  name.df$x = name.df$x
  
  p <- ggplot(plot.df2, mapping = aes_string(x = "x", y = "value")) +
    geom_line(mapping = aes_string(group = "variable", colour = "type"), size = 1.1,
              show.legend = show.legend) +
    #geom_text(name.df, mapping = aes_string(label = "variable", colour = "type",
    #                                       x = "x", y = "value"), show.legend = FALSE) +
    geom_label_repel(name.df, mapping = aes_string(label = "variable", colour = "type",
                                                   x = "x", y = "value"),max.overlaps = Inf, 
                     size = 5, show.legend = FALSE) +
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') +
    #annotate(geom="text", x = (min(plot.df2$x)+max(plot.df2$x))/2, y = max(plot.df2$value)*0.95,
    #         label=ind.name2, size = 6) + 
    scale_color_manual(
      name = "Effect Type",
      labels = c("Dep Effect", "All Effect", "Ind Effect"),
      values = darken(c("green", "blue", "red"))) +
    xlab("Habitat Index") + ylab("Niche Index") + xlab(NULL)+
    #ggtitle(ind.name2) + 
    theme_bw() + #ylab(NULL)+
    theme(plot.title = element_text(hjust = 0.5))+ theme(axis.text.y = element_text(hjust = 0))+ 
    #theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(-0.1,"cm"))+
    theme(plot.margin = unit(c(0,0,0,0),"lines"), 
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          axis.title = element_text(size = 10))
    
  
  return(p)
}

orig_df = df[which(rownames(df)==key_gene$names[select]),]


qdODE_plot_8 <- function(result,label = 10, j,show.legend = TRUE, nrow = NULL, ncol = NULL){
  get_a_plot <- function(j){

    p = list()
    for (i in c(6,4,7,2,5,3,1,8)) {
      p[[i]] = qdODE_plot_base(result[[i]], show.legend = FALSE)
      ra = ggplot_build(p[[i]])$layout$panel_params[[1]]$x.range
      p[[i]] = p[[i]]+scale_x_continuous(limits = c(ra[1],ra[2]),labels = math_format(expr = 10^.x))+
        theme(
          axis.title.x = element_text(size = 12),
          axis.text = element_text(size = 18),
          axis.title.y = element_text(size = 12))
        
      
    }
    
    for (i in c(1,2,5,6)) {
      p[[i]] = p[[i]]+ 
        annotate("rect",xmin=-Inf, xmax=Inf, ymin=-Inf,ymax=Inf, alpha=0.25, fill="#FF9BD2") 
      
    }
    
    for (i in c(3,4,7,8)) {
      p[[i]] = p[[i]]+ 
        annotate("rect",xmin=-Inf, xmax=Inf, ymin=-Inf,ymax=Inf, alpha=0.25, fill='#AAD7D9') 
      
    }
    
    y_limits = Reduce(rbind,lapply(1:8,function(c)ggplot_build(p[[c]])$layout$panel_params[[1]]$y.range))
    y_limits_lower = min(y_limits[,1])
    y_limits_upper = max(y_limits[,2]) 
    x_limits = Reduce(rbind,lapply(1:8,function(c)ggplot_build(p[[c]])$layout$panel_params[[1]]$x.range))
    x_limits_lower = min(x_limits[,1])
    x_limits_upper = max(x_limits[,2]) 
    
    for (i in 1:7) {
      p[[i]] = p[[i]] + scale_y_continuous(limits = c(y_limits_lower, y_limits_upper),
                                           labels = math_format(expr = 10^.x)) + 
        ylab(NULL) + scale_x_continuous(limits = c(x_limits_lower, x_limits_upper),
                                        labels = math_format(expr = 10^.x)) 
    }
    p[[8]] = p[[8]] +  scale_y_continuous(name = NULL, 
                                          limits = c(y_limits_lower, y_limits_upper),
                                          labels = math_format(expr = 10^.x)) + 
      scale_x_continuous(limits = c(x_limits_lower, x_limits_upper),
                         labels = math_format(expr = 10^.x)) 
    
    for (i in 2:8) {
      p[[i]] = p[[i]] + theme(axis.text.y = element_blank(), axis.ticks.length.y = unit(-0.1,"cm"))
    }
    
    

    
    pp = wrap_plots(p,ncol = 8)
    return(pp)
  }
  get_a_plot()
  
  pp = get_a_plot(i)

  #pp = wrap_plots(pp,nrow = 5)
  
  y_lab <- ggplot() + annotate(geom = "text", size=7,x = 1, y = 1,
                               label = "Gene Expression", angle = 90) +
    coord_cartesian(clip = "off") + theme_void()
  pp = (y_lab | wrap_plots(pp, nrow = 1, ncol = ncol)) +
    plot_annotation(caption = "Expression Index",
                    theme = theme(plot.caption = element_text(size = 18,hjust=.5))) +
    plot_layout(widths = c(.05, 1),guides = 'collect')
  
  return(pp)
}


Module_netM36 = lapply(1:8,function(c) lapply(ODE_submodule_all[[20]][[c]]$ode_result, network_conversion))

pdf("Fig8_A.pdf",width= 25,height= 30)
layout(matrix(1:8, 2, 4, byrow = TRUE))
lapply(1:8, function(c) network_plot(Module_netM36[[c]], title = NULL))
dev.off()


pdf("Fig8_B.pdf",width=15,height=10)
#vp1 <- viewport(height=unit(0.3, "npc"), width=unit(0.85, "npc"), y=0.15, x=0.5)
layout(matrix(1:8, 2, 4, byrow = TRUE))

network_plot1(Module_net[[6]], title = "SSS")
network_plot1(Module_net[[4]], title = "TTS")
network_plot1(Module_net[[7]], title = "STS")
network_plot1(Module_net[[2]], title = "TSS")

network_plot1(Module_net[[5]], title = "SSR")
network_plot1(Module_net[[3]], title = "TTR")
network_plot1(Module_net[[1]], title = "STR")
network_plot1(Module_net[[8]], title = "TSR")
#print(p1, vp=vp1)
dev.off()



pp = qdODE_plot_8(result = gene)
ggsave("Fig8C.pdf",pp,width = 40, height = 4)


