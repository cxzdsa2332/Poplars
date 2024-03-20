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
df <- read.csv("Code/graft_expression_matrix.csv",
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
  
  X = X[X!=0]
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



fit_all <- function(d1,d2,n1,n2){
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
  
  
  core.number <- 12
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {require(parallel)})
  clusterExport(cl, c("power_equation_fit", "power_fit", "df_new"), envir = environment())
  #all_model = parLapply(cl = cl, 1:length(new_df), function(c) power_equation_fit(new_df))
  all_model = parLapply(cl = cl,n1:n2, function(c) power_equation_fit(df_new[[c]]))
  stopCluster(cl)
  all_model
}

#lapply(1:10, function(c) power_equation_fit(df_new[[c]]))
result_fit_1 = list(fit_SS = fit_all(SSS,SSR,1,5),
                    fit_ST = fit_all(STS,STR,1,5),
                    fit_TS = fit_all(TSS,TSR,1,5),
                    fit_TT = fit_all(TTS,TTR,1,5))



get_interaction <- function(data, col, reduction = FALSE ){
  
  data <- t(data); name <- colnames(data)
  y = as.matrix(data[,col])
  x = as.matrix(data[,-col])
  
  #lm(y~x)
  n <- ncol(x)
  return_obj = list(ind.name = name[col],
                    dep.name = name[-col],
                    coefficient = as.numeric(coef(lm(y~x))[2]))
  
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
                                        parameters = qdODE.est$par))
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


qdODE_parallel <- function(result, reduction = FALSE, thread = 2, maxit = 1e3){
  data = result$original_data
  relationship = lapply(1:nrow(data),function(c)get_interaction(data, c, reduction = reduction))
  
  result = lapply(1:nrow(data),function(c) qdODE_all(result = result,
                                                     relationship = relationship,
                                                     i = c,
                                                     maxit = maxit
  ))
  names(result) = rownames(data)
  names(relationship) = rownames(data)
  return_obj <- list(ode_result = result,
                     relationship = relationship)
  return(return_obj)
}

qdODE_parallel2 <- function(result, thread = 2, maxit = 1e3){
  
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {require(orthopolynom)})
  clusterEvalQ(cl, {require(deSolve)})
  clusterExport(cl, c("qdODEmod", "qdODE_ls", "qdODE_fit", "qdODE_all","get_legendre_matrix",
                      "get_legendre_par","legendre_fit","qdODE_parallel",
                      "get_interaction"), envir=environment())
  result = parLapply(1:length(result),function(c) qdODE_parallel(result[[c]]), cl = cl)
  stopCluster(cl)
  
  return(result)
}


ode_all = list(qdODE_parallel2(result_fit_1$fit_SS, thread = 5),
               qdODE_parallel2(result_fit_1$fit_ST, thread = 5),
               qdODE_parallel2(result_fit_1$fit_TS, thread = 5),
               qdODE_parallel2(result_fit_1$fit_TT, thread = 5))



check <- function(j,s){
  effect_scion = mean(ode_all[[j]][[s]]$ode_result$Rootstock$fit$Scion)
  effect_rootstock = mean(ode_all[[j]][[s]]$ode_result$Scion$fit$Rootstock)
  
  
  threthod = 0.01
  
  if (abs(effect_scion)<=threthod) {
    effect_scion = 0
  }
  
  if (abs(effect_rootstock)<=threthod) {
    effect_rootstock = 0
  }
  
  
  
  if ( effect_scion > 0 & effect_rootstock > 0){
    type = "Dove-Dove"
  } else if (effect_scion == 0 & effect_rootstock > 0) {
    type = "Zero-Dove"
  } else if (effect_scion < 0 & effect_rootstock > 0) {
    type = "Hawk-Dove"
  } else if (effect_scion > 0 & effect_rootstock == 0) {
    type = "Dove-Zero"
  } else if (effect_scion == 0 & effect_rootstock == 0) {
    type = "Zero-Zero"
  } else if (effect_scion < 0 & effect_rootstock == 0) {
    type = "Hawk-Zero"
  } else if (effect_scion > 0 & effect_rootstock < 0) {
    type = "Dove-Hawk"
  } else if (effect_scion == 0 & effect_rootstock < 0) {
    type = "Zero-Hawk"
  } else if (effect_scion < 0 & effect_rootstock < 0) {
    type = "Hawk-Hawk"
  } 
  return(type)
}


check_all <- function(n){
  SS = data.frame(table(sapply(1:n,function(x)check(1,x))))
  ST = data.frame(table(sapply(1:n,function(x)check(2,x))))
  TS = data.frame(table(sapply(1:n,function(x)check(3,x))))
  TT = data.frame(table(sapply(1:n,function(x)check(4,x))))
  
  all_type = c("Dove-Dove","Zero-Dove","Hawk-Dove",
               "Dove-Zero","Zero-Zero","Hawk-Zero",
               "Dove-Hawk","Zero-Hawk","Hawk-Hawk")
  df_reference = data.frame(Var1 = all_type)
  df_reference = merge(df_reference,SS,by = "Var1",all.x = T)
  df_reference = merge(df_reference,ST,by = "Var1",all.x = T)
  df_reference = merge(df_reference,TS,by = "Var1",all.x = T)
  df_reference = merge(df_reference,TT,by = "Var1",all.x = T)
  colnames(df_reference) = c("type", "S/S","S/T","T/S","T/T")
  df_reference[is.na(df_reference)] = 0
  return(df_reference)
}

df_reference = check_all(nrow(df))



df$Stype = c("Dove","Dove","Dove","Hawk","Hawk","Hawk","Zero","Zero","Zero")
df$Rtype = c("Dove","Hawk","Zero","Dove","Hawk","Zero","Dove","Hawk","Zero")

df2 = melt(df)
df2$value = df2$value/sum(df2$value)*100


cols <- c("S.T" = '#FFB534', "T.S" = '#65B741', "S.S" = "#FF004D", "T.T" = "#3652AD")

p = ggplot(df2, mapping = aes(x = variable, y = value, fill = variable))+
  geom_bar(stat = 'identity') +
  facet_wrap(~factor(type, c("Dove-Dove","Zero-Dove","Hawk-Dove",
                             "Dove-Zero","Zero-Zero","Hawk-Zero",
                             "Dove-Hawk","Zero-Hawk","Hawk-Hawk")))+
  theme_bw()+
  theme(axis.title=element_text(size=15),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10,hjust = 0),
        panel.spacing = unit(0.0, "lines"),
        plot.margin = unit(c(1,1,1,1), "lines"),
        strip.background = element_blank(),
        plot.background = element_blank(),
        strip.text = element_blank())+
  #scale_y_continuous(expand = expansion(mult = 0, add = 0))+
  ylab("Hawk  Zero  Dove")+xlab("Dove  Zero  Hawk")+
  scale_fill_manual(values = cols)

