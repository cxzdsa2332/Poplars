rm(list=ls())
#load(file = "SAD8_31_40.Rdata")
library(mvtnorm)
library(parallel)
library(nlme)
library(drc)
library(aomisc)
library(patchwork)
library(ggplot2)
library(reshape2)
library(cowplot)
library(scales)

df = read.csv(file = "Code/graft_expression_matrix.csv", header = T, row.names = NULL)
df = df[,-c(2:5)] #remove control

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
df = Reduce(cbind, df_all2)


power_equation <- function(x, power_par){ t(sapply(1:nrow(power_par),
                                                   function(c) power_par[c,1]*x^power_par[c,2] ) )}
logsumexp <- function(v) {
  vm = max(v)
  log(sum(exp(v-vm))) + vm
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

power_equation_all <- function(x,y, maxit=1e2){
  result <- power_equation_base(x,y)
  iter <- 1
  while( is.null(result) && iter <= maxit) {
    iter <- iter + 1
    try(result <- power_equation_base(x,y))
  }
  return(result)
}

get_SAD1_covmatrix <- function(par,n){
  phi <- par[1]; gamma <- par[2];
  sigma <- array(dim=c(n,n))
  #formula 1, diag element
  diag(sigma) <- sapply(1:n, function(c)(1-phi^(2*c))/(1-phi^2) )
  #formula 2, non-diag element
  sigma[lower.tri(sigma)] <- do.call(c,lapply(1:(n-1),function(c)phi^seq(1:(n-c))*diag(sigma)[c]))
  sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(t(sigma))]
  return(gamma^2*sigma)
}

power_fit <- function(times,y){
  x <- as.numeric(times)
  y <- as.numeric(y)
  if (all(y==0)) {
    result = NULL
  } else{
    min_value = min(y[y!=0])
    
    lmFit <- lm( log( y + runif(1, min = 0, max = min_value))  ~ log(x))
    coefs <- coef(lmFit)
    a <- exp(coefs[1])
    b <- coefs[2]
    
    model <- try(nls(y~a*x^b,start = list(a = a, b = b),
                     control = nls.control(maxiter = 1e4, minFactor = 1e-200)))
    result = model
    if( 'try-error' %in% class(model)) {
      result = NULL
    }
  }  
  return(result)
}


power_fit2 <- function(times,y){
  tmp <- data.frame(x=as.numeric(times),y=as.numeric(y))  
  model <- try(drm(y ~ x, fct = DRC.powerCurve(),
                   data = tmp))
  x <- as.numeric(times)
  y <- as.numeric(y)
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a =0.3, b = 0.1),
                     control = nls.control(maxiter = 100000,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a =3, b = 0.1),
                     control = nls.control(maxiter = 100000,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a =0.1, b = 0.1),
                     control = nls.control(maxiter = 100000,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a = -0.1, b = 0.1),
                     control = nls.control(maxiter = 100000,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a = 0.8, b = -0.1),
                     control = nls.control(maxiter = 100000,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    model <- try(nls(y~a*x^b,start = list(a =-0.1, b = -0.1),
                     control = nls.control(maxiter = 1000000,minFactor = 1e-200)))
  }
  if ('try-error' %in% class(model)) {
    result <- c(NA,NA)
  }
  else{
    result <- coef(model)}
  return(result)
}

#EM
Q_function <- function(par, prob_log, omega_log, X, k, Time){
  n = dim(X)[1]; d = dim(X)[2];
  
  par.mu <- matrix(par[-c(1:16)],nrow = k,ncol = 2*8 )
  par.mu = list(par.mu[,1:2],par.mu[,3:4],par.mu[,5:6],par.mu[,7:8],par.mu[,9:10],
                par.mu[,11:12],par.mu[,13:14],par.mu[,15:16])
  
  
  par.cov = list(par[1:2],par[3:4],par[5:6],par[7:8],par[9:10],
                 par[11:12],par[13:14],par[15:16])
  cov = lapply(1:8,function(c) get_SAD1_covmatrix(par.cov[[c]], n = 8))
  mu <- lapply(1:8,function(c)power_equation(Time[[c]], par.mu[[c]]))
  
  mvn_log1 <- sapply(1:k, function(c) dmvnorm(X[,1:8],   mu[[1]][c,], cov[[1]], log = T))
  mvn_log2 <- sapply(1:k, function(c) dmvnorm(X[,9:16],  mu[[2]][c,], cov[[2]], log = T))
  mvn_log3 <- sapply(1:k, function(c) dmvnorm(X[,17:24], mu[[3]][c,], cov[[3]], log = T))
  mvn_log4 <- sapply(1:k, function(c) dmvnorm(X[,25:32], mu[[4]][c,], cov[[4]], log = T))
  mvn_log5 <- sapply(1:k, function(c) dmvnorm(X[,33:40], mu[[5]][c,], cov[[5]], log = T))
  mvn_log6 <- sapply(1:k, function(c) dmvnorm(X[,41:48], mu[[6]][c,], cov[[6]], log = T))
  mvn_log7 <- sapply(1:k, function(c) dmvnorm(X[,49:56], mu[[7]][c,], cov[[7]], log = T))
  mvn_log8 <- sapply(1:k, function(c) dmvnorm(X[,57:64], mu[[8]][c,], cov[[8]], log = T))
  
  mvn.log = mvn_log1 + mvn_log2 + mvn_log3 + mvn_log4 + mvn_log5 + mvn_log6 + mvn_log7 + mvn_log8
  
  tmp = sweep(mvn.log, 2, FUN = "+", STATS = prob_log) - omega_log
  Q = -sum(tmp*exp(omega_log))
  return(Q)
}

get_par_int8 <- function(X, k, times,Time){
  n = dim(X)[1]; d = dim(X)[2];
  
  cov.int = c(c(1e-5,mean(diag(cov(X[,1:8 ])))),
              c(1e-5,mean(diag(cov(X[,9:16])))),
              c(1e-5,mean(diag(cov(X[,17:24])))),
              c(1e-5,mean(diag(cov(X[,25:32])))),
              c(1e-5,mean(diag(cov(X[,33:40])))),
              c(1e-5,mean(diag(cov(X[,41:48])))),
              c(1e-5,mean(diag(cov(X[,49:56])))),
              c(1e-5,mean(diag(cov(X[,57:64])))))
  
  #cov.int = c(c(0.008,0.1),
  #            c(0.008,0.1),
  #            c(0.008,0.1),
  #            c(0.008,0.1),
  #            c(0.008,0.1),
  #            c(0.008,0.1),
  #            c(0.008,0.1),
  #            c(0.008,0.1))
  
  #init.cluster <- Mclust(X, k)
  #init.cluster$centers = t(init.cluster$parameters$mean)
  init.cluster <- kmeans(X,centers = k)
  prob <- table(init.cluster$cluster)/n
  
  
  fit1 <- lapply(1:k,function(c) power_fit2(Time[[1]],  init.cluster$centers[c,1:8]))
  fit2 <- lapply(1:k,function(c) power_fit2(Time[[2]], init.cluster$centers[c,9:16]))
  fit3 <- lapply(1:k,function(c) power_fit2(Time[[3]],init.cluster$centers[c,17:24]))
  fit4 <- lapply(1:k,function(c) power_fit2(Time[[4]],init.cluster$centers[c,25:32]))
  fit5 <- lapply(1:k,function(c) power_fit2(Time[[5]],init.cluster$centers[c,33:40]))
  fit6 <- lapply(1:k,function(c) power_fit2(Time[[6]],init.cluster$centers[c,41:48]))
  fit7 <- lapply(1:k,function(c) power_fit2(Time[[7]],init.cluster$centers[c,49:56]))
  fit8 <- lapply(1:k,function(c) power_fit2(Time[[8]],init.cluster$centers[c,57:64]))
  fit_all = list(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8)
  
  mu.par.int <- cbind(matrix(Reduce(c,fit1), k, 2, byrow = T),
                      matrix(Reduce(c,fit2), k, 2, byrow = T),
                      matrix(Reduce(c,fit3), k, 2, byrow = T),
                      matrix(Reduce(c,fit4), k, 2, byrow = T),
                      matrix(Reduce(c,fit5), k, 2, byrow = T),
                      matrix(Reduce(c,fit6), k, 2, byrow = T),
                      matrix(Reduce(c,fit7), k, 2, byrow = T),
                      matrix(Reduce(c,fit8), k, 2, byrow = T))
  
  par.mu <- mu.par.int
  par.mu = list(par.mu[,1:2],par.mu[,3:4],par.mu[,5:6],par.mu[,7:8],par.mu[,9:10],
                par.mu[,11:12],par.mu[,13:14],par.mu[,15:16])
  
  mu <- lapply(1:8,function(c)power_equation(Time[[c]], par.mu[[c]]))
  
  
  check.negative = function(fit, mu,i){
    realmu = init.cluster$centers[,(i*8-8+1):(i*8)]
    no = which( sapply(fit, "[[", 1) <0)
    no2 = which( abs(sapply(1:k, function(c) crossprod(mu[[i]][c,],realmu[c,]))) >=50)
    no = union(no,no2)
    if ( length(no) ==0 ) {
      fit = fit
    } else{
      tmp = lapply(no,function(c) power_fit(times[(i*8-8+1):(i*8)],  init.cluster$centers[c,(i*8-8+1):(i*8)]))
      fit[no] = lapply(tmp, coef)
    }
    return(fit)
  }
  
  fit_all2 = lapply(1:8, function(c) check.negative(fit_all[[c]],mu, c))
  
  
  mu.par.int2 <- cbind(matrix(Reduce(c,fit_all2[[1]]), k, 2, byrow = T),
                       matrix(Reduce(c,fit_all2[[2]]), k, 2, byrow = T),
                       matrix(Reduce(c,fit_all2[[3]]), k, 2, byrow = T),
                       matrix(Reduce(c,fit_all2[[4]]), k, 2, byrow = T),
                       matrix(Reduce(c,fit_all2[[5]]), k, 2, byrow = T),
                       matrix(Reduce(c,fit_all2[[6]]), k, 2, byrow = T),
                       matrix(Reduce(c,fit_all2[[7]]), k, 2, byrow = T),
                       matrix(Reduce(c,fit_all2[[8]]), k, 2, byrow = T))
  #fit1 <- lapply(1:k,function(c) power_fit(times[1:8],  init.cluster$centers[c,1:8]))
  #fit2 <- lapply(1:k,function(c) power_fit(times[9:16], init.cluster$centers[c,9:16]))
  #fit3 <- lapply(1:k,function(c) power_fit(times[17:24],init.cluster$centers[c,17:24]))
  #fit4 <- lapply(1:k,function(c) power_fit(times[25:32],init.cluster$centers[c,25:32]))
  #fit5 <- lapply(1:k,function(c) power_fit(times[33:40],init.cluster$centers[c,33:40]))
  #fit6 <- lapply(1:k,function(c) power_fit(times[41:48],init.cluster$centers[c,41:48]))
  #fit7 <- lapply(1:k,function(c) power_fit(times[49:56],init.cluster$centers[c,49:56]))
  #it8 <- lapply(1:k,function(c) power_fit(times[57:64],init.cluster$centers[c,57:64]))
  
  
  
  #mu.par.int <- cbind(t(sapply(fit1, coef)),
  #                    t(sapply(fit2, coef)),
  #                    t(sapply(fit3, coef)),
  #                    t(sapply(fit4, coef)),
  #                    t(sapply(fit5, coef)),
  #                    t(sapply(fit6, coef)),
  #                    t(sapply(fit7, coef)),
  #                    t(sapply(fit8, coef)))
  
  
  par.mu <- mu.par.int2
  par.mu = list(par.mu[,1:2],par.mu[,3:4],par.mu[,5:6],par.mu[,7:8],par.mu[,9:10],
                par.mu[,11:12],par.mu[,13:14],par.mu[,15:16])
  
  mu2 <- lapply(1:8,function(c)power_equation(Time[[c]], par.mu[[c]]))
  
  
  return_obj <- list(initial_cov_params = cov.int,
                     initial_mu_params = mu.par.int2,
                     initial_probibality = prob)
  return(return_obj)
}

#initial.pars = lapply(2:40,function(c)get_par_int8(X,c,times,Time))
#save.image(file = "SAD_initial_par.Rdata")
#data = df;k=10
fun_clu8 <- function(data, k, trans = log10, inv.cov = NULL,
                     initial.pars = NULL, iter.max = 1e2, parscale = 1e-20){
  df1 = data[,1:8]
  df1 = df1[,order(colSums(df1))]
  
  df2 = data[,9:16]
  df2 = df2[,order(colSums(df2))]
  
  df3 = data[,17:24]
  df3 = df3[,order(colSums(df3))]
  
  df4 = data[,25:32]
  df4 = df4[,order(colSums(df4))]
  
  df5 = data[,33:40]
  df5 = df5[,order(colSums(df5))]
  
  df6 = data[,41:48]
  df6 = df6[,order(colSums(df6))]
  
  df7 = data[,49:56]
  df7 = df7[,order(colSums(df7))]
  
  df8 = data[,57:64]
  df8 = df8[,order(colSums(df8))]
  
  #attribute
  data = as.matrix(data); n = dim(data)[1]; d = dim(data)[2]
  eplison = 1; iter = 0;
  
  times1 = as.numeric(trans(colSums(df1)+1))
  times2 = as.numeric(trans(colSums(df2)+1))
  times3 = as.numeric(trans(colSums(df3)+1))
  times4 = as.numeric(trans(colSums(df4)+1))
  times5 = as.numeric(trans(colSums(df5)+1))
  times6 = as.numeric(trans(colSums(df6)+1))
  times7 = as.numeric(trans(colSums(df7)+1))
  times8 = as.numeric(trans(colSums(df8)+1))
  
  X1 = trans(df1+1)
  X2 = trans(df2+1)
  X3 = trans(df3+1)
  X4 = trans(df4+1)
  X5 = trans(df5+1)
  X6 = trans(df6+1)
  X7 = trans(df7+1)
  X8 = trans(df8+1)
  
  times = c(times1,times2,times3,times4,times5,times6,times7,times8)
  Time = list(times1,times2,times3,times4,times5,times6,times7,times8)
  X = cbind(X1,X2,X3,X4,X5,X6,X7,X8)
  # initial pars
  if (is.null(initial.pars)) {
    initial.pars = get_par_int8(X, k, times,Time)
  }
  par.int <- c(initial.pars$initial_cov_params, initial.pars$initial_mu_params)
  prob_log <- log(initial.pars$initial_probibality)
  
  parscale = par.int[2]*parscale
  
  while( abs(eplison) > 1e-3 && iter <= iter.max ){
    #E step
    par.mu <- matrix(par.int[-c(1:16)],nrow = k,ncol = 2*8 )
    par.mu = list(par.mu[,1:2],par.mu[,3:4],par.mu[,5:6],par.mu[,7:8],par.mu[,9:10],
                  par.mu[,11:12],par.mu[,13:14],par.mu[,15:16])
    
    
    par.cov = list(par.int[1:2],par.int[3:4],par.int[5:6],par.int[7:8],par.int[9:10],
                   par.int[11:12],par.int[13:14],par.int[15:16])
    cov = lapply(1:8,function(c) get_SAD1_covmatrix(par.cov[[c]], n = 8))
    mu <- lapply(1:8,function(c)power_equation(Time[[c]], par.mu[[c]]))
    
    
    #par.cov[[c]
    mvn_log1 <- sapply(1:k, function(c) dmvnorm(X[,1:8],   mu[[1]][c,], cov[[1]], log = T))
    mvn_log2 <- sapply(1:k, function(c) dmvnorm(X[,9:16],  mu[[2]][c,], cov[[2]], log = T))
    mvn_log3 <- sapply(1:k, function(c) dmvnorm(X[,17:24], mu[[3]][c,], cov[[3]], log = T))
    mvn_log4 <- sapply(1:k, function(c) dmvnorm(X[,25:32], mu[[4]][c,], cov[[4]], log = T))
    mvn_log5 <- sapply(1:k, function(c) dmvnorm(X[,33:40], mu[[5]][c,], cov[[5]], log = T))
    mvn_log6 <- sapply(1:k, function(c) dmvnorm(X[,41:48], mu[[6]][c,], cov[[6]], log = T))
    mvn_log7 <- sapply(1:k, function(c) dmvnorm(X[,49:56], mu[[7]][c,], cov[[7]], log = T))
    mvn_log8 <- sapply(1:k, function(c) dmvnorm(X[,57:64], mu[[8]][c,], cov[[8]], log = T))
    
    mvn = mvn_log1 + mvn_log2 + mvn_log3 + mvn_log4 + mvn_log5 + mvn_log6 + mvn_log7 + mvn_log8
    mvn = sweep(mvn, 2, FUN = '+', STATS =  prob_log )
    
    
    #logsumexp(mvn[1,])
    omega_log = t(sapply(1:n, function(c) mvn[c,] - logsumexp(mvn[c,]) ))
    omega = exp(omega_log)
    
    LL.mem <- Q_function(par = par.int, prob_log = prob_log, omega_log, X, k, Time)
    
    #M step
    prob_exp = apply(omega_log, 2, logsumexp)
    prob_log = prob_exp - log(n)
    
    Q.maximization <- try(optim(par = par.int, Q_function,
                                prob_log = prob_log,
                                omega_log = omega_log,
                                X = X,
                                k = k,
                                Time = Time,
                                #method = "BFGS",
                                #lower = c(0,0,rep(-Inf, 2*k)),
                                #upper = Inf,
                                control = list(trace = TRUE,maxit = 500,
                                               parscale = c(rep(parscale,16),rep(1,2*k*8))
                                )))
    if ('try-error' %in% class(Q.maximization))
      break
    par.hat <- Q.maximization$par
    par.int = par.hat
    LL.next <- Q_function(par = par.hat, prob_log = prob_log, omega_log, X, k, Time)
    eplison <-  LL.next - LL.mem
    LL.mem <- LL.next
    iter = iter + 1
    
    cat("\n", "iter =", iter, "\n", "Log-Likelihood = ", LL.next, "\n")
  }
  AIC = 2*(LL.next) + 2*(length(par.hat)+k-1)
  BIC = 2*(LL.next) + log(n)*(length(par.hat)+k-1)
  
  omega = exp(omega_log)
  X.clustered <- data.frame(X, apply(omega,1,which.max), check.names = F)
  
  return_obj <- list(cluster_number = k,
                     Log_likelihodd = LL.mem,
                     AIC = AIC,
                     BIC = BIC,
                     cov_par = par.cov,
                     mu_par = par.mu,
                     probibality = exp(prob_log),
                     omega = omega,
                     cluster = X.clustered,
                     cluster2 = data.frame(data, apply(omega,1,which.max), check.names = F),
                     Time = times,
                     mu = mu, 
                     original_data = data)
  return(return_obj)
}


fun_clu_parallel <- function(data, Time = NULL, trans = log10, start, end, iter.max = 100, thread = 2){
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {library(mvtnorm)})
  clusterExport(cl,c(c("fun_clu8","initial.pars","Q_function","logsumexp","power_equation_base",
                       "power_equation_all","power_equation","get_SAD1_covmatrix"),ls()),
                envir=environment())
  result <- parLapply(cl=cl, start:end, function(c) fun_clu8(data = data,
                                                             trans = trans,
                                                             initial.pars = initial.pars[[c-1]],
                                                             k = c,
                                                             iter.max = iter.max))
  stopCluster(cl)
  
  return(result)
}

#data = df

#res = fun_clu8(data = df, k = 10, iter.max = 1,parscale = 1e-20)

result_31_40 = fun_clu_parallel(df,start = 2,end = 50,iter.max = 5,thread = 10)




result = result_31_40[[9]]

#fun_clu_BIC(result)




#result = result[[9]]
#best.k=10
#show.legend = T

cluster.result = result
times = result$Time
Time = list(times[1:8], times[9:16], times[17:24], times[25:32], 
            times[33:40], times[41:48], times[49:56], times[57:64])

no = as.numeric(names(table(result$cluster$`apply(omega, 1, which.max)`)))
mu_par = lapply(1:8,function(c) result$mu_par[[c]][no,])





ref2 <- data.frame(type = 1:8, 
                   Group = c("STR","TSS","TTR","TTS","SSR","SSS","STS","TSR"))

df_all = list(cluster.result$cluster[,c(1:8,65)],
              cluster.result$cluster[,c(9:16,65)],
              cluster.result$cluster[,c(17:24,65)],
              cluster.result$cluster[,c(25:32,65)],
              cluster.result$cluster[,c(33:40,65)],
              cluster.result$cluster[,c(41:48,65)],
              cluster.result$cluster[,c(49:56,65)],
              cluster.result$cluster[,c(57:64,65)])


get_plot_df <- function(Time0,mu_par0,df0,degree=0.1,i){
  times = Time0
  times_new = seq(min(times),max(times),length = 30)
  
  par.mu = mu_par0
  k = 36
  alpha = as.numeric(table(cluster.result$cluster$`apply(omega, 1, which.max)`))
  
  
  plot.df = df0

 

  plot.df$`apply(omega, 1, which.max)`[plot.df$`apply(omega, 1, which.max)`>=7] = 
    plot.df$`apply(omega, 1, which.max)`[plot.df$`apply(omega, 1, which.max)`>=7]-1
  

  plot.df$`apply(omega, 1, which.max)`[plot.df$`apply(omega, 1, which.max)`>=14] = 
    plot.df$`apply(omega, 1, which.max)`[plot.df$`apply(omega, 1, which.max)`>=14]-1
  

  
  
  plot.df$`apply(omega, 1, which.max)`[plot.df$`apply(omega, 1, which.max)`>=30] = 
    plot.df$`apply(omega, 1, which.max)`[plot.df$`apply(omega, 1, which.max)`>=30]-1
  
 # table(plot.df$`apply(omega, 1, which.max)`)

  plot.df$`apply(omega, 1, which.max)`
  
  X = plot.df[,-ncol(plot.df)]
  plot.df$name = rownames(plot.df)
  colnames(plot.df) = c(times,"cluster","name")
  
  plot.df = melt(plot.df, id.vars = c('cluster',"name"))
  colnames(plot.df) = c("cluster","name", "x","y")
  plot.df$x = as.numeric(as.character(plot.df$x))
  
  df.alpha = data.frame(cluster = 1:k, alpha = min(alpha)/alpha*degree)
  plot.df = merge(plot.df, df.alpha, by = "cluster")
  plot.df$type = i
  plot.df$name = paste0(plot.df$name,"_",plot.df$type)
  plot.df
}

plot.df_all = lapply(1:8,function(c) get_plot_df(Time[[c]], mu_par[[c]], df_all[[c]],i=c))
plot.df = Reduce(rbind, plot.df_all)


get_mu.fit <- function(Time0, mu_par0, i){
  times_new = seq(min(Time0), max(Time0), length = 100)
  mu.fit = power_equation(times_new, mu_par0)
  colnames(mu.fit) = times_new
  mu.fit = melt(as.matrix(mu.fit))
  colnames(mu.fit) = c("cluster","x","y")
  mu.fit$x = as.numeric(as.character(mu.fit$x))
  mu.fit$type = i
  mu.fit
}
mu.fit_all = lapply(1:8,function(c) get_mu.fit(Time[[c]], mu_par[[c]], i = c))
mu.fit = Reduce(rbind, mu.fit_all)
mu.fit = merge(mu.fit,ref2,by = "type",all.x = T)
mu.fit$linetype = substr(mu.fit$Group,3,3)
mu.fit$linetype[which(mu.fit$linetype=='R')]=1
mu.fit$linetype[which(mu.fit$linetype=='S')]=2


alpha = as.numeric(table(cluster.result$cluster$`apply(omega, 1, which.max)`))
name.df = data.frame(label = paste0("M",1:36," (",alpha ,")"),
                     x = mean(range(times)), y = max(cluster.result$cluster[,-65])*0.9, 
                     cluster = 1:36)

plot.df = merge(plot.df,ref2,by = "type",all.x = T)
#range(plot.df_all[[1]]$x)
#range(mu.fit_all[[6]]$x)


inter = c(which(plot.df$Group=="SSS"),
          which(plot.df$Group=="SSR"),
          which(plot.df$Group=="TTS"),
          which(plot.df$Group=="TTR"))
          
intra = c(which(plot.df$Group=="STS"),
          which(plot.df$Group=="STR"),
          which(plot.df$Group=="TSS"),
          which(plot.df$Group=="TSR"))

d1 = plot.df[inter,]
d2 = plot.df[intra,]


mu1 = c(which(mu.fit$Group=="SSS"),
        which(mu.fit$Group=="SSR"),
        which(mu.fit$Group=="TTS"),
        which(mu.fit$Group=="TTR"))
mu2 = c(which(mu.fit$Group=="STS"),
        which(mu.fit$Group=="STR"),
        which(mu.fit$Group=="TSS"),
        which(mu.fit$Group=="TSR"))

m1 = mu.fit[mu1,]
m2 = mu.fit[mu2,]


#i=1



cols <- c("STS" = '#FFB534', "TSS" = '#65B741', "SSS" = "#FF004D", "TTS" = "#3652AD",
          "STR" = '#FFB534', "TSR" = '#65B741', "SSR" = "#FF004D", "TTR" = "#3652AD")

#i=1
get_p1 <- function(i){
  e1 = d1[which(d1$cluster==i),]
  e2 = d2[which(d2$cluster==i),]
  
  mm1 = m1[which(m1$cluster==i),]
  mm2 = m2[which(m2$cluster==i),]
  
  p1 = ggplot() + geom_line(e1,
                            mapping = aes_string(x = "x", y = "y",
                                                 colour = factor(e1$Group),
                                                 group = "name"),
                            alpha = e1$alpha,
                            show.legend = F) +
    geom_line(mm1, mapping = aes_string(x = "x", y = "y", colour = factor(mm1$Group)),
              linetype = as.numeric(mm1$linetype),
              size=1.25, show.legend = F)+ 
    theme(axis.title=element_text(size=18)) +
    theme_bw() + 
    guides(alpha = "none") +
    scale_colour_manual(values = cols)+
    scale_y_continuous(limits = c(0,5),breaks = c(0:4),
                       labels = math_format(expr = 10^.x))+
    
    theme(axis.text.x = element_blank(), 
          axis.ticks.length.x = unit(-0.1,"cm"))+
    ylab(NULL)+xlab(NULL)+theme(plot.margin = unit(c(0,0,0,0),"lines"))
  
  
  p2 = ggplot() + geom_line(e2,
                            mapping = aes_string(x = "x", y = "y",
                                                 colour = factor(e2$Group),
                                                 group = "name"),
                            alpha = e2$alpha,
                            show.legend = F) +
    geom_line(mm2, mapping = aes_string(x = "x", y = "y",
                                        colour = factor(mm2$Group)),
              
              linetype = as.numeric(mm2$linetype),
              size=1.25, show.legend = F)+ 
    theme(axis.title=element_text(size=18)) +
    theme_bw() + 
    guides(alpha = "none") +
    scale_colour_manual(values = cols)+
    scale_y_continuous(limits = c(0,5),breaks = c(0:4))+
    
    theme(axis.text.x = element_blank(), 
          axis.ticks.length.x = unit(-0.1,"cm"),
          axis.text.y = element_blank(), 
          axis.ticks.length.y = unit(-0.1,"cm")) +
    ylab(NULL)+xlab(NULL)+
    theme(plot.margin = unit(c(0,0,0,0),"lines"))+
    ggtitle(name.df$label[i])
  p = p1+p2
  p
}

get_p2 <- function(i){
  e1 = d1[which(d1$cluster==i),]
  e2 = d2[which(d2$cluster==i),]
  
  mm1 = m1[which(m1$cluster==i),]
  mm2 = m2[which(m2$cluster==i),]
  
  p1 = ggplot() + geom_line(e1,
                            mapping = aes_string(x = "x", y = "y",
                                                 colour = factor(e1$Group),
                                                 group = "name"),
                            alpha = e1$alpha,
                            show.legend = F) +
    geom_line(mm1, mapping = aes_string(x = "x", y = "y", colour = factor(mm1$Group)),
              linetype = as.numeric(mm1$linetype),
              size=1.25, show.legend = F)+ 
    theme(axis.title=element_text(size=18)) +
    theme_bw() + 
    guides(alpha = "none") +
    scale_colour_manual(values = cols)+
    scale_y_continuous(limits = c(0,5),breaks = c(0:4))+
    
    theme(axis.text.x = element_blank(), 
          axis.ticks.length.x = unit(-0.1,"cm"),
          axis.text.y = element_blank(), 
          axis.ticks.length.y = unit(-0.1,"cm")) +
    ylab(NULL)+xlab(NULL)+theme(plot.margin = unit(c(0,0,0,0),"lines"))
  
  
  p2 = ggplot() + geom_line(e2,
                            mapping = aes_string(x = "x", y = "y",
                                                 colour = factor(e2$Group),
                                                 group = "name"),
                            alpha = e2$alpha,
                            show.legend = F) +
    geom_line(mm2, mapping = aes_string(x = "x", y = "y",
                                        colour = factor(mm2$Group)),
              
              linetype = as.numeric(mm2$linetype),
              size=1.25, show.legend = F)+ 
    theme(axis.title=element_text(size=18)) +
    theme_bw() + 
    guides(alpha = "none") +
    scale_colour_manual(values = cols)+
    scale_y_continuous(limits = c(0,5),breaks = c(0:4))+
    
    theme(axis.text.x = element_blank(), 
          axis.ticks.length.x = unit(-0.1,"cm"),
          axis.text.y = element_blank(), 
          axis.ticks.length.y = unit(-0.1,"cm")) +
    ylab(NULL)+xlab(NULL)+
    theme(plot.margin = unit(c(0,0,0,0),"lines"))+
    ggtitle(name.df$label[i])
  p = p1+p2
  p
}
  
get_p3 <- function(i){
  e1 = d1[which(d1$cluster==i),]
  e2 = d2[which(d2$cluster==i),]
  
  mm1 = m1[which(m1$cluster==i),]
  mm2 = m2[which(m2$cluster==i),]
  
  p1 = ggplot() + geom_line(e1,
                            mapping = aes_string(x = "x", y = "y",
                                                 colour = factor(e1$Group),
                                                 group = "name"),
                            alpha = e1$alpha,
                            show.legend = F) +
    geom_line(mm1, mapping = aes_string(x = "x", y = "y", colour = factor(mm1$Group)),
              linetype = as.numeric(mm1$linetype),
              size=1.25, show.legend = F)+ 
    theme(axis.title=element_text(size=18)) +
    theme_bw() + 
    guides(alpha = "none") +
    scale_colour_manual(values = cols)+
    scale_y_continuous(limits = c(0,5),breaks = c(0:4),labels = math_format(expr = 10^.x))+
    scale_x_continuous(labels = math_format(expr = 10^.x))+
    ylab(NULL)+xlab(NULL)+theme(plot.margin = unit(c(0,0,0,0),"lines"))
  
  
  p2 = ggplot() + geom_line(e2,
                            mapping = aes_string(x = "x", y = "y",
                                                 colour = factor(e2$Group),
                                                 group = "name"),
                            alpha = e2$alpha,
                            show.legend = F) +
    geom_line(mm2, mapping = aes_string(x = "x", y = "y",
                                        colour = factor(mm2$Group)),
              
              linetype = as.numeric(mm2$linetype),
              size=1.25, show.legend = F)+ 
    theme(axis.title=element_text(size=18)) +
    theme_bw() + 
    guides(alpha = "none") +
    scale_colour_manual(values = cols)+
    scale_y_continuous(limits = c(0,5),breaks = c(0:4))+
    scale_x_continuous(labels = math_format(expr = 10^.x))+
    
    theme(axis.text.y = element_blank(), 
          axis.ticks.length.y = unit(-0.1,"cm")) +
    ylab(NULL)+xlab(NULL)+
    theme(plot.margin = unit(c(0,0,0,0),"lines"))+
    ggtitle(name.df$label[i])
  p = p1+p2
  p
}
  
get_p4 <- function(i){
  e1 = d1[which(d1$cluster==i),]
  e2 = d2[which(d2$cluster==i),]
  
  mm1 = m1[which(m1$cluster==i),]
  mm2 = m2[which(m2$cluster==i),]
  
  p1 = ggplot() + geom_line(e1,
                            mapping = aes_string(x = "x", y = "y",
                                                 colour = factor(e1$Group),
                                                 group = "name"),
                            alpha = e1$alpha,
                            show.legend = F) +
    geom_line(mm1, mapping = aes_string(x = "x", y = "y", colour = factor(mm1$Group)),
              linetype = as.numeric(mm1$linetype),
              size=1.25, show.legend = F)+ 
    theme(axis.title=element_text(size=18)) +
    theme_bw() + 
    guides(alpha = "none") +
    scale_colour_manual(values = cols)+
    scale_y_continuous(limits = c(0,5),breaks = c(0:4))+
    
    theme(axis.text.x = element_blank(), 
          axis.ticks.length.x = unit(-0.1,"cm")) +
    scale_x_continuous(labels = math_format(expr = 10^.x))+
    ylab(NULL)+xlab(NULL)+theme(plot.margin = unit(c(0,0,0,0),"lines"))
  
  
  p2 = ggplot() + geom_line(e2,
                            mapping = aes_string(x = "x", y = "y",
                                                 colour = factor(e2$Group),
                                                 group = "name"),
                            alpha = e2$alpha,
                            show.legend = F) +
    geom_line(mm2, mapping = aes_string(x = "x", y = "y",
                                        colour = factor(mm2$Group)),
              
              linetype = as.numeric(mm2$linetype),
              size=1.25, show.legend = F)+ 
    theme(axis.title=element_text(size=18)) +
    theme_bw() + 
    guides(alpha = "none") +
    scale_colour_manual(values = cols)+
    scale_y_continuous(limits = c(0,5),breaks = c(0:4))+
    scale_x_continuous(labels = math_format(expr = 10^.x))+
    theme(axis.text.x = element_blank(), 
          axis.ticks.length.x = unit(-0.1,"cm")) +
    ylab(NULL)+xlab(NULL)+
    theme(plot.margin = unit(c(0,0,0,0),"lines"))+
    ggtitle(name.df$label[i])
  p = p1+p2
  p
}
  
  
a = 1:33

p1 = lapply(1:33,function(c)get_p3(c))
p1[[1]]


#p = list(p1[[1]],p2[[1]],p2[[2]]])
#p[[1]]
pp = cowplot::plot_grid(plotlist=p1,ncol = 3)
ggsave("Fig5.pdf",pp,width = 20 ,height = 30)


