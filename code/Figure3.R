rm(list=ls())
library(reshape2)
library(ggplot2)
library(nlme)
library(drc)
library(aomisc)
library(patchwork)
library(cowplot)
library(parallel)
library(ggrepel)
library(scales)
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
  x <- as.numeric(times)
  y <- as.numeric(y)
  
  y[y==0] = y +  min(y[y!=0])
  x[x==0] = x +  min(x[x!=0])
  lmFit <- lm( log( y )  ~ log(x))
  coefs <- coef(lmFit)
  a <- exp(coefs[1])
  b <- coefs[2]
  
  model <- try(nls(y~a*x^b,start = list(a = a, b = b),
                   control = nls.control(maxiter = 1e3, minFactor = 1e-200)))
  if( 'try-error' %in% class(model)) {
    tmp <- data.frame(x=as.numeric(times),y=as.numeric(y))  
    model <- try(drm(y ~ x, fct = DRC.powerCurve(),
                     data = tmp))
  } 
  if ('try-error' %in% class(model)) {
    result <- NA
  }
  else{
    result <- model}
  return(result)
}

df_R = df2[,grep("R",colnames(df2))]
df_S = df2[,-grep("R",colnames(df2))]


power_equation_fit <- function(data, n=30, trans = log10, thread = 2) {
  data = data[,order(colSums(data))]
  if ( is.null(trans)) {
    X = colSums(data)
    trans_data = data
  } else{
    X = trans(colSums(data)+1)
    trans_data = trans(data+1)
  }
  colnames(trans_data) = X
  
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {require(drc)})
  clusterEvalQ(cl, {require(aomisc)})
  clusterEvalQ(cl, {library(nlme)})
  clusterExport(cl, c("power_fit","trans_data", "X"), envir = environment())
  all_model = parLapply(cl = cl, 1:nrow(data), function(c) power_fit(trans_data[c,],X))
  stopCluster(cl)
  
  
  
  names(all_model) = rownames(data)
  no = which(sapply(all_model, length)>=1)
  all_model2 = all_model[no]
  data2 = data[no,]
  trans_data2 = trans_data[no,]
  
  new_x = seq(min(X), max(X), length = n)
  power_par = t(vapply(all_model2, coef, FUN.VALUE = numeric(2), USE.NAMES = TRUE))
  power_fit = t(vapply(all_model2, predict, newdata = data.frame(x=new_x),
                       FUN.VALUE = numeric(n), USE.NAMES = TRUE))
  
  colnames(power_fit) = new_x
  result = list(original_data = data2, trans_data = trans_data2,
                power_par = power_par, power_fit = power_fit,
                Time = X)
  return(result)
}

#power_equation_fit(df_new[[1]])


SSS = df_S[,grep(paste0("SSS"),colnames(df_S))]
STS = df_S[,grep(paste0("STS"),colnames(df_S))]
TSS = df_S[,grep(paste0("TSS"),colnames(df_S))]
TTS = df_S[,grep(paste0("TTS"),colnames(df_S))]

SSR = df_R[,grep(paste0("SSR"),colnames(df_R))]
STR = df_R[,grep(paste0("STR"),colnames(df_R))]
TSR = df_R[,grep(paste0("TSR"),colnames(df_R))]
TTR = df_R[,grep(paste0("TTR"),colnames(df_R))]



fit_all <- function(d1,d2,n){
  df_all = cbind(d1,d2)
  
  df_new = lapply(1:nrow(df_all),function(x) 
    data.frame(rbind(as.numeric(df_all[x,1:8]),as.numeric(df_all[x,9:16]))))
  f_rename = function(x){
    rownames(x)= c("Scion", "Rootstock")
    x
  }
  df_new = lapply(df_new, f_rename)
  
  
  names(df_new) = rownames(df_all)
  times = lapply(df_new,function(x) log10(colSums(x)+1))
  df_new2 = lapply(df_new,function(x)log10(x+1))
  
  
  core.number <- 6
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {require(parallel)})
  clusterExport(cl, c("power_equation_fit", "power_fit", "df_new"), envir = environment())
  #all_model = parLapply(cl = cl, 1:length(new_df), function(c) power_equation_fit(new_df))
  all_model = parLapply(cl = cl,1:n, function(c) power_equation_fit(df_new[[c]]))
  stopCluster(cl)
  all_model
}

#lapply(1:10, function(c) power_equation_fit(df_new[[c]]))
result_fit = list(fit_SS = fit_all(SSS,SSR,10),
                  fit_ST = fit_all(STS,STR,10),
                  fit_TS = fit_all(TSS,TSR,10),
                  fit_TT = fit_all(TTS,TTR,10))



get_interaction <- function(data, col, reduction = FALSE ){
  require(glmnet)
  
  data <- t(data); name <- colnames(data)
  y = as.matrix(data[,col])
  x = as.matrix(data[,-col])
  
  #lm(y~x)
  n <- ncol(x)
  
  
  ridge_cv <- try(cv.glmnet(x = x, y = y, type.measure = "mse", nfolds = 10, alpha = 0))
  if ('try-error' %in% class(ridge_cv)) {
    return_obj = list(ind.name = name[col],
                      dep.name = name[-col],
                      coefficient = as.numeric(coef(lm(y~x))[2]))
    
    
  } else{
    #plot(ridge_cv)
    best_ridge_coef <- abs(as.numeric(coef(ridge_cv, s = ridge_cv$lambda.min))[-1])
    
    fit <- cv.glmnet(x = x, y = y, alpha = 1, family = "gaussian", type.measure = "mse",
                     penalty.factor = 1/best_ridge_coef,
                     nfolds = 10, keep = TRUE
                     #, thresh=1e-30, maxit=1e6
    )
    lasso_coef <- coef(fit, s = fit$lambda.min)
    
    return_obj = list(ind.name = name[col],
                      dep.name = lasso_coef@Dimnames[[1]][lasso_coef@i + 1][-1],
                      coefficient = lasso_coef@x[-1])
    
    if ( length(return_obj$dep.name)==0 ) {
      lasso_coef <- coef(fit, s = fit$lambda.min)
      return_obj = list(ind.name = name[col],
                        dep.name = lasso_coef@Dimnames[[1]][lasso_coef@i + 1][-1],
                        coefficient = lasso_coef@x[-1])
    }
  }
  
  
  return(return_obj)
}

#get_interaction(data = r2$original_data,col=5)
darken <- function(color, factor=1.2){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

get_legendre_matrix <- function(x,legendre_order){
  legendre_coef <- legendre.polynomials(n = legendre_order, normalized=F)
  legendre_matrix <- as.matrix(as.data.frame(polynomial.values(
    polynomials = legendre_coef, x = scaleX(x, u = -1, v = 1))))
  colnames(legendre_matrix) <- paste0("legendre_",0:legendre_order)
  return(legendre_matrix[,2:(legendre_order+1)])
}

get_legendre_par <- function(y,legendre_order,x){
  #lm_method
  legendre_par <- as.numeric(coef(lm(y~get_legendre_matrix(x,legendre_order))))
  return(legendre_par)
}

legendre_fit <- function(par,x){
  legendre_order = length(par)
  fit <- sapply(1:length(par), function(c)
    par[c]*legendre.polynomials(n=legendre_order, normalized=F)[[c]])
  legendre_fit <- as.matrix(as.data.frame(polynomial.values(
    polynomials = fit, x = scaleX(x, u = -1, v = 1))))
  x_interpolation <- rowSums(legendre_fit)
  return(x_interpolation)
}

qdODEmod <- function(Time, State, Pars, power_par) {
  nn = length(Pars)
  ind_effect = paste0("alpha","*",names(State)[1])
  dep_effect = sapply(2:nn, function(c) paste0(paste0("beta",c-1),"*",names(State)[c]))
  dep_effect = paste0(dep_effect, collapse = "+")
  all_effect = paste0(ind_effect, "+", dep_effect)
  expr = parse(text = all_effect)
  
  with(as.list(c(State, Pars)), {
    dx = eval(expr)
    dy <- power_par[,1]*power_par[,2]*Time^(power_par[,2]-1)
    dind = alpha*x
    for(i in c(1:(nn-1))){
      tmp = paste0(paste0("beta",i),"*",paste0("y",i))
      expr2 = parse(text = tmp)
      assign(paste0("ddep",i),eval(expr2))
    }
    return(list(c(dx, dy, dind, mget(paste0("ddep",1:(nn-1))))))
  })
}

qdODE_ls <- function(pars, data, Time, power_par){
  n = length(pars)
  power_par = as.matrix(power_par)
  if (n==2) {
    Pars = c(alpha = pars[1], beta1 = pars[2:n])
    power_par = t(power_par)
    State = c(x=data[1,1],y1 = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  } else{
    Pars = c(alpha = pars[1], beta = pars[2:n])
    State = c(x=data[1,1],y = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  }
  out = as.data.frame(ode(func = qdODEmod, y = State, parms = Pars,
                          times = Time, power_par = power_par))
  X = as.numeric(data[1,])
  fit = as.numeric(out[,2])
  ind = as.numeric(out[,(n+2)])
  sse = sum(crossprod(X-fit),sum((ind[ind<0])^2))
  return(sse)
}

qdODE_fit <- function(pars, data, Time, power_par, LOP_order = 6, new_time = NULL, n_expand = 100){
  n = length(pars)
  if (n==2) {
    Pars = c(alpha = pars[1], beta1 = pars[2:n])
    power_par = t(power_par)
    State = c(x=data[1,1],y1 = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  } else{
    Pars = c(alpha = pars[1], beta = pars[2:n])
    State = c(x=data[1,1],y = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  }
  out = as.data.frame(ode(func = qdODEmod, y = State, parms = Pars,
                          times = Time, power_par = power_par))
  out2 = data.frame(x = out[,1], y = data[1,], y.fit = out[,2],
                    ind = out[,(n+2)], dep = out[,(n+3):(ncol(out))])
  colnames(out2)[4:ncol(out2)] = c(rownames(data)[1], rownames(data)[2:n])
  rownames(out2) = NULL
  
  all_LOP_par = sapply(2:ncol(out2),function(c)get_legendre_par(out2[,c], LOP_order, out2$x))
  
  if (is.null(new_time)) {
    time2 = seq(min(Time), max(Time), length = n_expand)
    out3 = apply(all_LOP_par, 2, legendre_fit, x = time2)
    out3 = cbind(time2, out3)
  } else{
    out3 = apply(all_LOP_par, 2, legendre_fit, x = new_time)
    out3 = cbind(new_time, out3)
  }
  colnames(out3) = colnames(out2)
  result = list(fit = out2,
                predict = data.frame(out3),
                LOP_par = all_LOP_par)
  return(result)
}

#pars = pars_int
#relationship = result2$relationship;result = result1[[2]];i=11;init_pars = 1;maxit = 1e5;pars = qdODE.est$par
qdODE_all <- function(result, relationship, i, init_pars = 1, LOP_order = 6, method = "ls",
                      new_time = NULL, n_expand = 100, maxit = 1e3){
  Time = as.numeric(colnames(result$power_fit))
  variable = c(relationship[[i]]$ind.name, relationship[[i]]$dep.name)
  data = result$power_fit[variable,]
  
  if (length(variable)<=1) {
    qdODE.est = NA
    result = NA
    return.obj <- append(result, list(ODE.value = NA,
                                      parameters = NA))
  } else{
    power_par = result$power_par[variable,][-1,]
    n = nrow(data)
    pars_int = c(init_pars,relationship[[i]]$coefficient)
    if (method == "ls") {
      qdODE.est <- optim(pars_int, qdODE_ls, data = data, Time = Time, power_par = power_par,
                         method = "L-BFGS-B",
                         lower = c(rep(-10,(length(pars_int)))),
                         upper = c(rep(10,(length(pars_int)))),
                         control = list(trace = TRUE, maxit = maxit))
      
      result <- qdODE_fit(pars = qdODE.est$par,
                          data = data,
                          power_par = power_par,
                          Time = Time)
      
      return.obj <- append(result, list(ODE.value = qdODE.est$value,
                                        parameters = qdODE.est$par,
                                        original_data = log10(result$original_data+1)))
    } else{
      qdODE.est <- optim(pars_int, qdODE_ls, data = data, Time = Time, power_par = power_par,
                         method = "L-BFGS-B",
                         lower = c(rep(-10,(length(pars_int)))),
                         #lower = c(0, rep(-10,(length(pars_int))-1)),
                         upper = c(rep(10,(length(pars_int)))),
                         control = list(trace = TRUE, maxit = maxit))
      
      result <- qdODE_fit(pars = qdODE.est$par,
                          data = data,
                          power_par = power_par,
                          Time = Time)
      return.obj <- append(result, list(ODE.value = qdODE.est$value,
                                        parameters = qdODE.est$par))
    }
  }
  return(return.obj)
}


#result = result_fit[[1]][[1]]
#result$original_data
qdODE_parallel <- function(result, reduction = FALSE, thread = 2, maxit = 1e3){
  data = result$original_data
  relationship = lapply(1:nrow(data),function(c)get_interaction(data, c, reduction = reduction))
  
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {require(orthopolynom)})
  clusterEvalQ(cl, {require(deSolve)})
  clusterExport(cl, c("qdODEmod", "qdODE_ls", "qdODE_fit", "qdODE_all","get_legendre_matrix",
                      "get_legendre_par","legendre_fit"), envir=environment())
  result = parLapply(1:nrow(data),function(c) qdODE_all(result = result,
                                                        relationship = relationship,
                                                        i = c,
                                                        maxit = maxit
  ), cl = cl)
  stopCluster(cl)
  names(result) = rownames(data)
  names(relationship) = rownames(data)
  return_obj <- list(ode_result = result,
                     relationship = relationship)
  return(return_obj)
}

#result=fit_SS[[1]];reduction = FALSE;thread = 2;maxit = 1e4



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

qdODE_plot_base <- function(result,original_data,time,s){
  result2 = qdODEplot_convert(result)
  plot.df2 = result2$plot.df2
  name.df = result2$name.df
  ind.name2 = result2$ind.name2
  name.df = name.df[-1,]
  name.df$x = name.df$x
  
  
  tmp = data.frame(x = as.numeric(time), y = as.numeric(log10(original_data[s,]+1)))
  
  if(s==1){
    shape=1
  } else{
    shape = 20
  }
  
  p <- ggplot(plot.df2, mapping = aes_string(x = "x", y = "value")) +
    geom_line(mapping = aes_string(group = "variable", colour = "type"), size = 1.1,
              show.legend = F) +
    geom_point(tmp, mapping = aes(x,y), shape = shape, alpha = 0.85,color = "#38419D") +
    #geom_text(name.df, mapping = aes_string(label = "variable", colour = "type",
    #                                       x = "x", y = "value"), show.legend = FALSE) +
    #geom_label_repel(name.df, mapping = aes_string(label = "variable", colour = "type",
    #                                               x = "x", y = "value"), show.legend = FALSE) +
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') +
    #annotate(geom="text", x = (min(plot.df2$x)+max(plot.df2$x))/2, y = max(plot.df2$value)*0.95,
    #         label=ind.name2, size = 6) + 
    scale_color_manual(
      #name = "Effect Type",
      #labels = c("Dep Effect", "All Effect", "Ind Effect"),
      values = c("#65B741", "#38419D", "#B80000")) +
    xlab("Habitat Index") + ylab("Niche Index") + 
    #ggtitle(ind.name2) + 
    theme_bw() + xlab(NULL)+ ylab(NULL)+ theme(legend.position = "none")+
    theme(plot.title = element_text(hjust = 0.5))+ theme(axis.text.y = element_text(hjust = 0))+ 
    #theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(-0.1,"cm"))+
    theme(plot.margin = unit(c(0,0,0,0),"lines"), 
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          axis.title = element_text(size = 10))#+ expand_limits(x=c(5.54,5.645))
  
  return(p)
}

#result = ode_all[[1]]$ode_result$Scion
#original_data=original_dataall[[1]]
#time=time_all[[1]]
#j=1

qdODE_plot_combined <- function(result,original_data,time,j){
  
  d1 = result$ode_result$Scion$fit
  n1 = ncol(d1)
  
  d2 = result$ode_result$Rootstock$fit
  n2 = ncol(d2)
  
  xm = mean(d1$Time)
  
  lower1 = min(d1[,-1])
  upper1 = max(d1[,-1])
  lower2 = min(d2[,-1])
  upper2 = max(d2[,-1])
  
  y_min = floor(min(lower1, lower2))
  y_max = ceiling(max(upper1, upper2))
  
  p1 = qdODE_plot_base(result$ode_result$Scion,original_data,time,s = 1) +
    scale_y_continuous(limits = c(y_min,y_max),labels = math_format(expr = 10^.x)) +
    theme(axis.text.y = element_text(hjust = 0))+
    scale_x_continuous(labels = math_format(expr = 10^.x))+ylab(NULL)+xlab(NULL)+
    theme(axis.text = element_text(size =10))
  
  p2 =  qdODE_plot_base(result$ode_result$Rootstock,original_data,time,s = 2) +
    scale_y_continuous(limits = c(y_min,y_max),sec.axis = sec_axis( trans=~.*10, name=rownames(df2)[j])) + 
    theme(axis.title.y.left = element_blank(),
          axis.text.y.left = element_blank(), 
          axis.text.y.right = element_blank(), 
          axis.ticks.y.right = element_blank(),
          axis.ticks.length.y = unit(-0.1,"cm"),
          axis.title.y = element_text(size = 15),
          axis.title.y.right = element_text(angle=-90, vjust = 0.5))+
    theme(axis.text = element_text(size =10))+
    theme(plot.title = element_text(hjust = 0.5,size = 18))+ 
    scale_x_continuous(labels = math_format(expr = 10^.x))+xlab(NULL)
  
  p = p1+p2
  
  p
}

plot_gene <- function(result_fit,j){
  
  ode_all = list(qdODE_parallel(result_fit$fit_SS[[j]]),
                 qdODE_parallel(result_fit$fit_ST[[j]]),
                 qdODE_parallel(result_fit$fit_TS[[j]]),
                 qdODE_parallel(result_fit$fit_TT[[j]]))
  
  original_dataall = list(result_fit$fit_SS[[j]]$original_data,
                          result_fit$fit_ST[[j]]$original_data,
                          result_fit$fit_TS[[j]]$original_data,
                          result_fit$fit_TT[[j]]$original_data)
  
  
  time_all = list(result_fit$fit_SS[[j]]$Time,
                  result_fit$fit_ST[[j]]$Time,
                  result_fit$fit_TS[[j]]$Time,
                  result_fit$fit_TT[[j]]$Time)
  
  p1 = qdODE_plot_combined(ode_all[[1]],original_dataall[[1]],time_all[[1]],j)
  p2 = qdODE_plot_combined(ode_all[[2]],original_dataall[[2]],time_all[[2]],j)
  p3 = qdODE_plot_combined(ode_all[[3]],original_dataall[[3]],time_all[[3]],j)
  p4 = qdODE_plot_combined(ode_all[[4]],original_dataall[[4]],time_all[[4]],j)
  
  pp = plot_grid(p1, p3, p2, p4)
  pp
  
}

p1 = plot_gene(result_fit, 1)
p2 = plot_gene(result_fit, 10)

ppp = p1/p2
ppp


ggsave("Fig3.pdf",width = 12, height = 9)
