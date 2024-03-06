rm(list=ls())
library(reshape2)
library(ggplot2)
library(nlme)
library(drc)
library(nlraa)
library(aomisc)
library(patchwork)
library(cowplot)
library(ggrepel)
library(scales)
df = read.csv(file = "rawdata/graft_expression_matrix.csv", header = T, row.names = 1)
df = df[,-c(1:4)]

df_ck = read.csv("Code/graft_expression_matrix.csv",
                 header = T,row.names = 1)[,1:4]
df_ck =df_ck[rownames(df),]



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

df_R = df[,grep("R",colnames(df))]
df_S = df[,-grep("R",colnames(df))]



d1 = data.frame(x = rep(c(5,6,8,12,18,26,36,48),4),
                y = log10(colSums(df_R)+1))
d1$type = c(rep("STR",8),rep("TTR",8),rep("SSR",8),rep("TSR",8))
d1$type3 = c(rep("S/T",8),rep("T/T",8),rep("S/S",8),rep("T/S",8))
d2 = data.frame(x = rep(c(5,6,8,12,18,26,36,48),4),
                y = log10(colSums(df_S)+1))
d2$type = c(rep("TSS",8),rep("TTS",8),rep("SSS",8),rep("STS",8))
d2$type3 = c(rep("T/S",8),rep("T/T",8),rep("S/S",8),rep("S/T",8))

d = rbind(d1,d2)
d$type2= c(rep("R",32),rep("S",32))






#i=1

get_p <- function(i){
  d_title = d[d$x==48,]
  d_title$x = d_title$x+0.1
  d_title = list(d_title[c(7,8),],
                 d_title[c(3,4),],
                 d_title[c(1,2),],
                 d_title[c(5,6),])
  
  cols <- c("S/T" = '#FFB534', "T/S" = '#65B741', "S/S" = "#FF004D", "T/T" = "#3652AD")
  
  ddd = list(rbind(d[d$type=="SSS",],d[d$type=="STS",]),
             rbind(d[d$type=="SSR",],d[d$type=="TSR",]),
             rbind(d[d$type=="STR",],d[d$type=="TTR",]),
             rbind(d[d$type=="TSS",],d[d$type=="TTS",]))
  df = ddd[[i]]
  d_title = d_title[[i]]
  
  
  
  new_d = data.frame(x=c(0,0,0,0),y=log10(colSums(df_ck+1)),
                     type = colnames(df_ck),
                     type3 = colnames(df_ck),
                     type2 = c("R","S","R","S"))
  
  
  
  tmp1 = new_d[4,]
  #tmp1 = rbind(tmp1,tmp1)
  tmp1$type = c("SSS")
  tmp1$type3 = c("S/S")
  
  tmp3 = new_d[3,]
  #tmp3 = rbind(tmp3)
  tmp3$type = c("SSR")
  tmp3$type3 = c("S/S")
  
  tmp2 = new_d[1,]
  #tmp2 = rbind(tmp2,tmp2)
  tmp2$type = c("TTR")
  tmp2$type3 = c("T/T")
  
  tmp4 = new_d[2,]
  #tmp4 = rbind(tmp4,tmp4)
  tmp4$type = c("TTS")
  tmp4$type3 = c("T/T")
  
  
  tmp_all = list(tmp1,tmp3,tmp2,tmp4)
  
  
  df2 = rbind(df,tmp_all[[i]])
  
  ylabel2 = parse(text= paste(10,"^", c(5.54,5.58,5.62), sep="") )
  
  if (i==1) {
    p = ggplot()+
      geom_line(df2,mapping = aes(x,y,group = type,color = type3),
                size = 1.25,linetype = 2)+
      geom_label_repel(d_title,mapping = aes(x,y,label = type3,color = type3))+
      theme(legend.position = "none")+theme_bw()+
      scale_y_continuous(limits = c(5.54,5.63),breaks = c(5.54,5.58,5.62),
                         labels = ylabel2)+
      scale_x_continuous(breaks = c(5,10,20,30,40,50))+
      xlab(NULL)+ylab(NULL)+
      scale_colour_manual(values = cols)+
      scale_fill_manual(values =cols)
  } 
  if (i==2) {
    p = ggplot()+
      geom_line(df2,mapping = aes(x,y,group = type,color = type3),
                size = 1.25,linetype = 1)+
      geom_label_repel(d_title,mapping = aes(x,y,label = type3,color = type3))+
      theme(legend.position = "none")+theme_bw()+
      scale_y_continuous(limits = c(5.54,5.63),breaks = c(5.54,5.58,5.62),
                         labels = ylabel2)+
      scale_x_continuous(breaks = c(5,10,20,30,40,50))+
      xlab(NULL)+ylab(NULL)+
      scale_colour_manual(values = cols)+
      scale_fill_manual(values =cols)
  }
  if (i==3) {
    p = ggplot()+
      geom_line(df2,mapping = aes(x,y,group = type,color = type3),
                size = 1.25,linetype = 1)+
      geom_label_repel(d_title,mapping = aes(x,y,label = type3,color = type3))+
      theme(legend.position = "none")+theme_bw()+
      scale_y_continuous(limits = c(5.54,5.63),breaks = c(5.54,5.58,5.62),
                         labels = ylabel2)+
      scale_x_continuous(breaks = c(5,10,20,30,40,50))+
      xlab(NULL)+ylab(NULL)+
      scale_colour_manual(values = cols)+
      scale_fill_manual(values =cols)
  }
  if (i==4) {
    p = ggplot()+
      geom_line(df2,mapping = aes(x,y,group = type,color = type3),
                size = 1.25,linetype = 2)+
      geom_label_repel(d_title,mapping = aes(x,y,label = type3,color = type3))+
      theme(legend.position = "none")+theme_bw()+
      scale_y_continuous(limits = c(5.54,5.63),breaks = c(5.54,5.58,5.62),
                         labels = ylabel2)+
      scale_x_continuous(breaks = c(5,10,20,30,40,50))+
      xlab(NULL)+ylab(NULL)+
      scale_colour_manual(values = cols)+
      scale_fill_manual(values =cols)
  }
  
  
  #p = p+xlab(NULL)+ylab(NULL)+
  #scale_y_continuous(labels = ylabel2)#+
  #scale_x_reverse()+coord_flip()
  p
}

pp = lapply(1:4,function(x)get_p(x))

p = wrap_plots(pp,ncol = 1)
p

ggsave("fig1C.pdf",width = 6, height = 5)





power_equation_plot <- function(n = 8){
  set.seed(2022)
  
  df = new_df[,-1]
  no = sample(1:nrow(df),n)
  #no = c(1,2,3,4,5,6,7,8,9,10)
  data1 = df
  
  colnames(data1) = log10(colSums(data1)+1)
  
  df_original_1 =  as.matrix(log10(data1[no,]+1))
  Times1 = log10(colSums(data1)+1)
  
  mod1 = lapply(1:n,function(c) power_equation_base(x = Times1[1:8], y = df_original_1[c,1:8]))
  mod2 = lapply(1:n,function(c) power_equation_base(x = Times1[9:16], y = df_original_1[c,9:16]))
  mod3 = lapply(1:n,function(c) power_equation_base(x = Times1[17:24], y = df_original_1[c,17:24]))
  mod4 = lapply(1:n,function(c) power_equation_base(x = Times1[25:32], y = df_original_1[c,25:32]))
  mod5 = lapply(1:n,function(c) power_equation_base(x = Times1[33:40], y = df_original_1[c,33:40]))
  mod6 = lapply(1:n,function(c) power_equation_base(x = Times1[41:48], y = df_original_1[c,41:48]))
  mod7 = lapply(1:n,function(c) power_equation_base(x = Times1[49:56], y = df_original_1[c,49:56]))
  mod8 = lapply(1:n,function(c) power_equation_base(x = Times1[57:64], y = df_original_1[c,57:64]))
  
  
  no2 = Reduce(intersect, list(which(sapply(mod1,is.null)==FALSE),
                               which(sapply(mod2,is.null)==FALSE),
                               which(sapply(mod3,is.null)==FALSE),
                               which(sapply(mod4,is.null)==FALSE),
                               which(sapply(mod5,is.null)==FALSE),
                               which(sapply(mod6,is.null)==FALSE),
                               which(sapply(mod7,is.null)==FALSE),
                               which(sapply(mod8,is.null)==FALSE)))
  
  no2 = c(56,137,168,193)
  
  df_original_1 =  df_original_1[no2,]
  
  d1 = data.frame(df_original_1)
  dd = data.frame(melt(df_original_1))
  dd2 = dd[order(dd$Var1,decreasing = F),]
  rownames(dd2) = NULL
  
  #dd2$group = c(rep("STR",32),rep("TTS",32),rep("TTR",32),rep("TTS",32),
  #             rep("SSR",32),rep("SSS",32),rep("STS",32),rep("TSR",32))
  
  
  
  
  dat_fit1 = lapply(mod1[no2], predict2_nls, interval = "confidence")
  dat_fit2 = lapply(mod2[no2], predict2_nls, interval = "confidence")
  dat_fit3 = lapply(mod3[no2], predict2_nls, interval = "confidence")
  dat_fit4 = lapply(mod4[no2], predict2_nls, interval = "confidence")
  dat_fit5 = lapply(mod5[no2], predict2_nls, interval = "confidence")
  dat_fit6 = lapply(mod6[no2], predict2_nls, interval = "confidence")
  dat_fit7 = lapply(mod7[no2], predict2_nls, interval = "confidence")
  dat_fit8 = lapply(mod8[no2], predict2_nls, interval = "confidence")
  
  
  for (i in 1:4 ){
    #tmp1 = data.frame(row.names(df_original_1)[i], Times1[1:8] )
    dat_fit1[[i]] = cbind( data.frame(row.names(df_original_1)[i], Times1[1:8] ),dat_fit1[[i]])
    dat_fit2[[i]] = cbind( data.frame(row.names(df_original_1)[i], Times1[9:16] ),dat_fit2[[i]])
    dat_fit3[[i]] = cbind( data.frame(row.names(df_original_1)[i], Times1[17:24] ),dat_fit3[[i]])
    dat_fit4[[i]] = cbind( data.frame(row.names(df_original_1)[i], Times1[25:32] ),dat_fit4[[i]])
    dat_fit5[[i]] = cbind( data.frame(row.names(df_original_1)[i], Times1[33:40] ),dat_fit5[[i]])
    dat_fit6[[i]] = cbind( data.frame(row.names(df_original_1)[i], Times1[41:48] ),dat_fit6[[i]])
    dat_fit7[[i]] = cbind( data.frame(row.names(df_original_1)[i], Times1[49:56] ),dat_fit7[[i]])
    dat_fit8[[i]] = cbind( data.frame(row.names(df_original_1)[i], Times1[57:64] ),dat_fit8[[i]])
    
    colnames(dat_fit1[[i]])[1:3] = c("Var1","Var2","value")
    colnames(dat_fit2[[i]])[1:3] = c("Var1","Var2","value")
    colnames(dat_fit3[[i]])[1:3] = c("Var1","Var2","value")
    colnames(dat_fit4[[i]])[1:3] = c("Var1","Var2","value")
    colnames(dat_fit5[[i]])[1:3] = c("Var1","Var2","value")
    colnames(dat_fit6[[i]])[1:3] = c("Var1","Var2","value")
    colnames(dat_fit7[[i]])[1:3] = c("Var1","Var2","value")
    colnames(dat_fit8[[i]])[1:3] = c("Var1","Var2","value")
    
  }
  
  
  
  df.fit1 = do.call(rbind,dat_fit1)
  df.fit2 = do.call(rbind,dat_fit2)
  df.fit3 = do.call(rbind,dat_fit3)
  df.fit4 = do.call(rbind,dat_fit4)
  df.fit5 = do.call(rbind,dat_fit5)
  df.fit6 = do.call(rbind,dat_fit6)
  df.fit7 = do.call(rbind,dat_fit7)
  df.fit8 = do.call(rbind,dat_fit8)
  
  df.fit1$group = "STR"
  df.fit2$group = "TSS"
  df.fit3$group = "TTR"
  df.fit4$group = "TTS"
  df.fit5$group = "SSR"
  df.fit6$group = "SSS"
  df.fit7$group = "STS"
  df.fit8$group = "TSR"
  
  #display.brewer.all()
  #df.text = data.frame(Var1 = rownames(data1[no2,][c(5,11,13,20,23,33,35,45,58, 
  #                                                   61,68,74),]), x = 5.6, y = 2.4)
  
  df_original_2 = rbind(df.fit1,df.fit2,df.fit3,df.fit4,
                        df.fit5,df.fit6,df.fit7,df.fit8)
  
  df_original_3 = df_original_2[order(df_original_2$Var1,decreasing = T),]
  df_original_3 = cbind(df_original_3,c(dd2$value[1:64],dd2$value[193:256],dd2$value[129:192],dd2$value[65:128]))
  colnames(df_original_3)[8] = c("original")
  
  #col = brewer.pal(n = 7, name = "Dark2")
  p3 <- ggplot() +
    geom_point(df_original_3, mapping = aes(x = Var2, y = original, group = group,color = group),
               size = 1.2, show.legend = F, alpha = 0.5) +
    geom_line(df_original_3, mapping = aes(x = Var2, y = value, group = group,color = group),
              size = 1.2, show.legend = F, alpha = 0.5) +
    geom_ribbon(df_original_3, mapping = aes(x = Var2, y = value, ymin = Q2.5, ymax = Q97.5,fill = group), 
                alpha = 0.2) +
    facet_grid(vars(group),vars(Var1))+
    
    
    xlab("Expression Index") + ylab("Individual Expression") + theme(axis.title=element_text(size=18)) +
    scale_x_continuous(labels = math_format(expr = 10^.x)) +
    scale_y_continuous(labels = math_format(expr = 10^.x)) +
    theme_bw() +
    theme(axis.title=element_text(size=15),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10,hjust = 0),
          #panel.spacing = unit(0.0, "lines"),
          plot.margin = unit(c(1,1,1,1), "lines"),
          strip.background = element_blank(),
          plot.background = element_blank(),
          #strip.text = element_blank()
    )
  return(p)
}

power_equation_plot()



#dff = log10(colSums(df)+1)
dff = colSums(df)+1
bardata <- data.frame(SSR = dff[grep("SSR",names(dff))],SSS = dff[grep("SSS",names(dff))],
                      TTR = dff[grep("TTR",names(dff))],TTS = dff[grep("TTS",names(dff))],
                      STR = dff[grep("STR",names(dff))],STS = dff[grep("STS",names(dff))],
                      TSR = dff[grep("TSR",names(dff))],TSS = dff[grep("TSS",names(dff))])



bardata2 = data.frame(id = substr(colnames(bardata),1,2), 
                      pos = substr(colnames(bardata),3,3), 
                      mean = apply(bardata,2,mean),
                      sd = apply(bardata,2,sd))


p1<-ggplot(bardata2, aes(x=id, y=mean, fill=pos)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean +sd),
                position=position_dodge(.8), width=.2) +
  theme_bw()+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim = c(0, 4.5e05))+
  theme(axis.text.x = element_text(size = 14, color = "black"))+##设置x轴字体大小
  theme(axis.text.y = element_text(size = 14, color = "black"))+##设置y轴字体大小
  theme(title=element_text(size=13))+#设置标题字体大小
  theme_bw()+xlab("Graft Combinations")+ylab("Niche Index")


p1
###############


power_equation <- function(x, power_par){ t(sapply(1:nrow(power_par),
                                                   function(c) power_par[c,1]*x^power_par[c,2] ) )}

power_equation_fit <- function(data, n = 30){
  data = as.matrix(data)
  EI = log10(colSums(data)+1)
  data = log10(data+1)
  model = lapply(1:nrow(data),function(c) power_equation_base(EI, data[c,]))
  confidence = lapply(model, predict2_nls, interval = "confidence")
  
  new_EI = seq(min(EI),max(EI),length=n)
  
  df_fit = sapply(1:2,function(c) predict(model[[c]], newdata = data.frame(x = new_EI)))
  rownames(df_fit) = new_EI
  colnames(df_fit) = c("R","S")
  dd = melt(df_fit)
  dd$type = substr(c(rep(rownames(data)[1],30),rep(rownames(data)[2],30)),1,2)
  
  confidence[[1]]$type = rownames(data)[1]
  confidence[[2]]$type = rownames(data)[2]
  confidence2 = rbind(confidence[[1]],confidence[[2]])
  
  orig_dd = data.frame(x = EI, y1 = data[1,], y2 = data[2,])
  colnames(orig_dd) = c("x",substr(rownames(data),3,3))
  orig_dd = melt(orig_dd, id.var = "x")
  tmp = cbind(orig_dd, confidence2)
  tmp[,8] = substr(tmp[,8],1,2)
  return(list(dd,tmp))
}

power_equation_base <- function(x, y){
  x <- as.numeric(x)
  y <- as.numeric(y)
  min_value = min(y[y!=0])
  
  lmFit <- lm( log( y + runif(1, min = 0, max = min_value))  ~ log(x))
  coefs <- coef(lmFit)
  a <- exp(coefs[1])
  b <- coefs[2]
  
  model <- try(nls(y~a*x^b,start = list(a = a, b = b),
                   control = nls.control(maxiter = 1e3, minFactor = 1e-200)))
  if( 'try-error' %in% class(model)) {
    result = NULL
  }
  else{
    result = model
  }
  return(result)
}


ds = data.frame(t(bardata))
colnames(ds) = NULL
ds_list = list(ds[1:2,],ds[3:4,],ds[5:6,],ds[7:8,])
all = lapply(1:4,function(c)power_equation_fit(ds_list[[c]]))


plt_df_orig = Reduce(rbind,lapply(all,'[[',2))

plt_df_fit = Reduce(rbind,lapply(all,'[[',1))
colnames(plt_df_fit)=c("x","variable","value","type")





y=colSums(df_ck+1)
yT = log10(y[1]+y[2])
yS = log10(y[3]+y[4])
y = log10(y)

new_d = data.frame(x=c(yT,yT,yS,yS),value=y,
                   variable = c("R","S","R","S"),
                   type = c("SS","SS","TT","TT"))





p2 = ggplot() + geom_point(plt_df_orig, mapping = aes(x=x,y=value,color=variable))+
  geom_ribbon(plt_df_orig, mapping = aes(x = x, y = value, ymin = Q2.5, ymax = Q97.5,fill = variable), 
              alpha = 0.2) +
  geom_point(new_d, mapping = aes(x=x,y=value,color=variable),shape = 2)+
  geom_line(plt_df_fit, mapping = aes(x=x,y=value,color=variable),linewidth=1.15)+
  facet_wrap(~type)+theme_bw()+xlab("Niche Index")+ylab("Organ Expression")+theme_bw() +
  scale_x_continuous(labels = math_format(expr = 10^.x)) +
  scale_y_continuous(labels = math_format(expr = 10^.x)) +
  theme(axis.title=element_text(size=15),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10,hjust = 0),
        #panel.spacing = unit(0.0, "lines"),
        plot.margin = unit(c(1,1,1,1), "lines"),
        strip.background = element_blank(),
        plot.background = element_blank(),
        #strip.text = element_blank()
  )
p2
#pp = ((p1 / p2) | p3 ) + plot_layout(guides = 'collect')+ plot_annotation(tag_levels = 'A')
ggsave("Fig1D.pdf",p2,width = 8, height = 5)
