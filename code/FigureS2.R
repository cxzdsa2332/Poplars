rm(list = ls())
library(igraph)
library(qgraph)
library(grid)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(reshape2)
load("Module_gene_ODE.Rdata") #result from ODE solving

darken <- function(color, factor=1.2){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

normalization <- function(x, z = 0.2){(x-min(x))/(max(x)-min(x))+z}

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
                                                   x = "x", y = "value"), show.legend = FALSE) +
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') +
    annotate(geom="text", x = (min(plot.df2$x)+max(plot.df2$x))/2, y = max(plot.df2$value)*0.95,
             label=ind.name2, size = 6) + 
    scale_color_manual(
      name = "Effect Type",
      labels = c("Dep Effect", "All Effect", "Ind Effect"),
      values = darken(c("green", "blue", "red"))) +
    xlab("Habitat Index") + ylab("Niche Index") + 
    #ggtitle(ind.name2) + 
    theme_bw() + xlab(NULL)+
    theme(plot.title = element_text(hjust = 0.5))+ theme(axis.text.y = element_text(hjust = 0))+ 
    #theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(-0.1,"cm"))+
    theme(plot.margin = unit(c(0,0,0,0),"lines"), 
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          axis.title = element_text(size = 10))#+ expand_limits(x=c(5.54,5.645))
  
  return(p)
}
#result= result3$ode_result$M1
#qdODE_plot_base(result3$ode_result$M1)
qdODE_plot_all <- function(result,label = 10, show.legend = TRUE, nrow = NULL, ncol = NULL){
  p = lapply(result$ode_result, qdODE_plot_base, label = label, show.legend = show.legend)
  #p = lapply(1:length(result$ode_result), function(c) qdODE_plot_base(result$ode_result[[c]], label = label))
  p = lapply(p, "+", xlab(NULL))
  p = lapply(p, "+", ylab(NULL))
  
  y_lab <- ggplot() + annotate(geom = "text", size=7,x = 1, y = 1,
                               label = "Module Expression", angle = 90) +
    coord_cartesian(clip = "off") + theme_void()
  pp = (y_lab | wrap_plots(p, nrow = nrow, ncol = ncol)) +
    plot_annotation(caption = "Expression Index",
                    theme = theme(plot.caption = element_text(size = 20,hjust=.5))) +
    plot_layout(widths = c(.05, 1),guides = 'collect')
  
  return(pp)
}


qdODE_plot_base(result = gene_M10$ode_result$MblContig100407)

qdODE_plot_all(result = ODE_M19)

pp = lapply(1:100, function(c) qdODE_plot_base(result = gene_M17$ode_result[[c]]))
pp = wrap_plots(pp)
ggsave(pp, filename = "1.pdf",width = 45, height = 45)



p1 = qdODE_plot_base(result = gene_M10$ode_result$MblContig100407)
p2 = qdODE_plot_base(result = gene_M10$ode_result[[10]])
p3 = qdODE_plot_base(result = ODE_M19$ode_result[[2]])
p4 = qdODE_plot_base(result = ODE_M19$ode_result[[4]])
pp = p1+p2+p3+p4


ggsave(pp, filename = "FigS2.pdf",width = 12, height = 7)
