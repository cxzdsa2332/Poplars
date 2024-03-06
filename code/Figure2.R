rm(list=ls())
library(reshape2)
library(ggplot2)
library(nlme)
library(drc)
library(aomisc)
library(patchwork)
library(cowplot)
df <- read.csv("rawdata/graft_expression_matrix.csv",
               header = T,row.names = NULL)

df = df[,-c(2:5)]
data_cleaning <- function(data, x = round(ncol(data)*0.3)){
  data = aggregate(data[,2:ncol(data)], by=list(data[,1]), FUN = 'sum')
  rownames(data) = data[,1]
  data = data[,-1]
  tmp = apply(data, 1, function(c) which( as.numeric(c) != 0) )
  keep_no = which(sapply(tmp, length) >= x)
  data2 = data[keep_no,]
  return(data2)
}



df_all = list(df[,c(1:9)],df[,c(1,10:17)],df[,c(1,18:25)],df[,c(1,26:33)],
              df[,c(1,34:41)],df[,c(1,42:49)],df[,c(1,50:57)],df[,c(1,58:65)])
df_all = lapply(df_all, data_cleaning,x=4)



common_gene = Reduce(intersect,lapply(df_all, rownames))

df_all2 = lapply(df_all,function(x) x[common_gene,])
df2 = Reduce(cbind, df_all2)


power_equation <- function(x, power_par){ t(sapply(1:nrow(power_par),
                                                   function(c) power_par[c,1]*x^power_par[c,2] ) )}


power_fit <- function(y,times){
  tmp <- data.frame(x=as.numeric(times),y=as.numeric(y))  
  model <- try(drm(y ~ x, fct = DRC.powerCurve(),
                   data = tmp))
  x <- as.numeric(times)
  y <- as.numeric(y)
  model <- try(nls(y~a*x^b,start = list(a =5, b = 0.5),
                   control = nls.control(maxiter = 1e5,minFactor = 1e-200)))
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a =0.3, b = 0.1),
                     control = nls.control(maxiter = 1e5,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a =3, b = 0.1),
                     control = nls.control(maxiter = 1e5,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a =0.1, b = 0.1),
                     control = nls.control(maxiter = 1e5,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a = -0.1, b = 0.1),
                     control = nls.control(maxiter = 1e5,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a = 0.8, b = -0.1),
                     control = nls.control(maxiter = 1e5,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a =-0.1, b = -0.1),
                     control = nls.control(maxiter = 1e50,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a =1e-100, b = 100),
                     control = nls.control(maxiter = 1e50,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a =1e100, b = -100),
                     control = nls.control(maxiter = 1e50,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    y = rep(mean(y),length(y))
    tmp <- data.frame(x=as.numeric(times),y=as.numeric(y))  
    model <- try(drm(y ~ x, fct = DRC.powerCurve(),
                     data = tmp))
  }
  if ('try-error' %in% class(model)) {
    result <- c(a=0,b=0)
  }
  else{
    result <- model}
  return(result)
}

df_R = df2[,grep("R",colnames(df2))]
df_S = df2[,-grep("R",colnames(df2))]


get_plot_data <- function(ddd,ID){
  d1 = ddd[,grep(paste0("SS",ID),colnames(ddd))]
  d2 = ddd[,grep(paste0("ST",ID),colnames(ddd))]
  d3 = ddd[,grep(paste0("TS",ID),colnames(ddd))]
  d4 = ddd[,grep(paste0("TT",ID),colnames(ddd))]
  set.seed(2023)
  n = sample(nrow(d1),8)
  
  dd = list(log10(d1[n,]+1),
            log10(d2[n,]+1),
            log10(d3[n,]+1),
            log10(d4[n,]+1))
  times = list(as.numeric(log10(colSums(d1)+1)), 
               as.numeric(log10(colSums(d2)+1)), 
               as.numeric(log10(colSums(d3)+1)), 
               as.numeric(log10(colSums(d4)+1)))
  
  m1 = lapply(1:length(n), function(c) power_fit(y = dd[[1]][c,], times = times[[1]]))
  m2 = lapply(1:length(n), function(c) power_fit(y = dd[[2]][c,], times = times[[2]]))
  m3 = lapply(1:length(n), function(c) power_fit(y = dd[[3]][c,], times = times[[3]]))
  m4 = lapply(1:length(n), function(c) power_fit(y = dd[[4]][c,], times = times[[4]]))
  
  
  new_times = lapply(1:4,function(c) seq(min(times[[c]]),max(times[[c]]),length=30))
  
  df_fit1 = sapply(1:length(n),function(c) predict(m1[[c]], newdata = data.frame(x = new_times[[1]])))
  colnames(df_fit1) = rownames(dd[[1]])
  df_fit2 = sapply(1:length(n),function(c) predict(m2[[c]], newdata = data.frame(x = new_times[[2]])))
  colnames(df_fit2) = rownames(dd[[2]])
  df_fit3 = sapply(1:length(n),function(c) predict(m3[[c]], newdata = data.frame(x = new_times[[3]])))
  colnames(df_fit3) = rownames(dd[[3]])
  df_fit4 = sapply(1:length(n),function(c) predict(m4[[c]], newdata = data.frame(x = new_times[[4]])))
  colnames(df_fit4) = rownames(dd[[4]])
  
  df_fitall = list(df_fit1,df_fit2,df_fit3,df_fit4)
  df_all = dd
  return(list(df_all,df_fitall,times,new_times))
}

result_rootstock = get_plot_data(df_R,ID = "R")
result_scion = get_plot_data(df_S,ID = "S")

result = result_rootstock
#i=1
get_base_plot <- function(i){
  data1 = cbind(as.numeric(result_rootstock[[1]][[1]][i,]),
                as.numeric(result_rootstock[[1]][[2]][i,]),
                as.numeric(result_rootstock[[1]][[3]][i,]),
                as.numeric(result_rootstock[[1]][[4]][i,]))
  
  
  colnames(data1) = c("SS","ST","TS","TT")
  
  data2 = cbind(result_rootstock[[2]][[1]][,i],
                result_rootstock[[2]][[2]][,i],
                result_rootstock[[2]][[3]][,i],
                result_rootstock[[2]][[4]][,i])
  colnames(data2) = c("SS","ST","TS","TT")
  data1 = melt(data1)
  data2 = melt(data2)
  data1$Var1 = unlist(result_rootstock[[3]])
  data2$Var1 = unlist(result_rootstock[[4]])
  
  
  data3 = cbind(as.numeric(result_scion[[1]][[1]][i,]),
                as.numeric(result_scion[[1]][[2]][i,]),
                as.numeric(result_scion[[1]][[3]][i,]),
                as.numeric(result_scion[[1]][[4]][i,]))
  
  
  colnames(data3) = c("SS","ST","TS","TT")
  
  data4 = cbind(result_scion[[2]][[1]][,i],
                result_scion[[2]][[2]][,i],
                result_scion[[2]][[3]][,i],
                result_scion[[2]][[4]][,i])
  colnames(data4) = c("SS","ST","TS","TT")
  data3 = melt(data3)
  data4 = melt(data4)
  data3$Var1 = unlist(result_scion[[3]])
  data4$Var1 = unlist(result_scion[[4]])
  
  data1$type = "Rootstock"
  data3$type = "Scion"
  data2$type = "Rootstock"
  data4$type = "Scion"
  
  d1 = rbind(data1,data3)
  d2 = rbind(data2,data4)
  
  
  cols <- c("ST" = '#FFB534', "TS" = '#65B741', "SS" = "#FF004D", "TT" = "#3652AD")
  cols2 <- c("ST" = '#AAD7D9', "TS" = '#AAD7D9', "SS" = "#FF9BD2", "TT" = "#FF9BD2")
  p = ggplot()+ 
    geom_rect(data = d1,aes(fill = Var2),xmin = -Inf,xmax = Inf,
              ymin = -Inf,ymax = Inf,alpha = 0.01) +
    geom_point(d1, mapping = aes(Var1, value,color = Var2,shape = type),alpha = 0.35)+
    geom_line(d2, mapping = aes(Var1, value,color = Var2,group = type,linetype=type),
              size = 0.85)+
    facet_wrap(~Var2)+
    theme_bw() +# scale_color_manual(values = c("Rootstock" = "#A94438", "Scion" = 'blue')) + 
    theme(axis.title=element_text(size=15),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10,hjust = 0),
          panel.spacing = unit(0.0, "lines"),
          plot.margin = unit(c(1,1,1,1), "lines"),
          strip.background = element_blank(),
          plot.background = element_blank(),
          strip.text = element_blank())+
    scale_color_manual(values = cols) + 
    scale_fill_manual(values =cols2)+
    scale_shape_manual(values = c(1,2))
  
  
  xlabel = ggplot_build(p)$layout$panel_params[[1]]$x.sec$breaks
  ylabel = ggplot_build(p)$layout$panel_params[[1]]$y.sec$breaks
  xlabel2 = parse(text= paste(10,"^", xlabel, sep="") )
  ylabel2 = parse(text= paste(10,"^", ylabel, sep="") )
  
  
  
  
  p = p+xlab(NULL)+ylab(NULL)+ggtitle(rownames(result[[1]][[1]])[i])+
    scale_x_continuous(labels = xlabel2) + 
    scale_y_continuous(labels = ylabel2)
  
  
  
  #p = ggplot()+ geom_point(data1, mapping = aes(Var1, value,color = Var2),shape = 1)+
  #  geom_line(data2, mapping = aes(Var1, value,color = Var2))+
  #  facet_wrap(~Var2)+
  #  theme_bw() + 
  #  theme(axis.title=element_text(size=15),
  #        axis.text.x = element_text(size=10),
  #        axis.text.y = element_text(size=10,hjust = 0),
  #        panel.spacing = unit(0.0, "lines"),
  #        plot.margin = unit(c(1,1,1,1), "lines"),
  #        strip.background = element_blank(),
  #        plot.background = element_blank(),
  #        strip.text = element_blank())
  #
  #xlabel = ggplot_build(p)$layout$panel_params[[1]]$x.sec$breaks
  #ylabel = ggplot_build(p)$layout$panel_params[[1]]$y.sec$breaks
  #xlabel2 = parse(text= paste(10,"^", xlabel, sep="") )
  #ylabel2 = parse(text= paste(10,"^", ylabel, sep="") )
  
  #p = p+scale_x_continuous(labels = xlabel2) + 
  #  scale_y_continuous(labels = ylabel2)+ 
  #  theme(panel.background = element_rect(fill = col))+
  #  ggtitle(rownames(result[[1]][[i]])[i])+
  # xlab(NULL)+ylab(NULL)
  
  #+
  #xlab(NULL)+ylab(NULL)+theme(panel.grid.major=element_line(colour="#DCF2F1"),
  #                            panel.grid.minor=element_line(colour="#DCF2F1"))
  
  p
}
pp = lapply(1:8,function(c) get_base_plot(c))

pp = wrap_plots(pp,ncol = 2, nrow = 4) + plot_layout(guides = 'collect')
ggsave("Fig2.pdf",pp,width = 10, height = 12)

