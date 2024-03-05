rm(list = ls())
library(ggplot2)
library(patchwork)
library(parallel)

nn = 30

#1. data generate by game model--------------------

real_par1 = c(0.01, 0.05, 0.1)
real_par2 = c(-0.15, 0.05, 0.25)
real_par3 = c(-0.55, -0.25, 0.35)


GLV <- function(t, x, parameters){
  with(as.list(c(x, parameters)), {
    #x[x < 10^-8] <- 0 # prevent numerical problems
    dxdt <- x * r + A %*% x
    list(dxdt)
  })
}

GLV3 <- function(t, x, parameters, dep_list){  
  with(as.list(c(x, parameters)), {
    n=3
    dxdt <- x[1:n] * r + A %*% x[1:n]
    
    for(i in c(1:n)){
      tmp = paste0(paste0("r","[",i,"]"),"*",paste0("x","[",i,"]"))
      expr2 = parse(text = tmp)
      assign(paste0("ind",i),eval(expr2))
    }
    tmp = list()
    for (i in 1:n) {
      tmp[[i]] = A[i,dep_list[[i]]]*x[dep_list[[i]]]
    }
    tmp = unlist(tmp)
    list(c(dxdt,mget(paste0("ind",1:n)), tmp))
  })
}

real_par = c(real_par1, real_par2, real_par3)
par = real_par

real_data <-  function(par){
  r = c(par[1],par[4],par[7])
  A = matrix(c(0, par[2], par[3] ,par[5], 0, par[6], par[8], par[9], 0), 
             nrow = 3, ncol = 3, byrow = T)
  x0 <- c(1, 1, 2)
  times = seq(4,7,length = nn)
  
  parameters <- list(r = r, A = A)
  # solve numerically
  out <- ode(y = x0, times = times, 
             func = GLV, parms = parameters)[,-1]
  
  dep_list = list(c(2,3), c(1,3), c(1,2))
  
  x00 = c(x0,x0,rep(0,6))
  out2 = data.frame(ode(y = x00, times = times, 
                        func = GLV3, parms = parameters,
                        dep_list = dep_list))
  
  out = data.frame(out)
  colnames(out) = c("x","y","z")
  colnames(out2) = c("times","x","y","z","x_ind","y_ind","z_ind", 
                     "y>-x","z>-x",
                     "x>-y","z>-y",
                     "x>-z","y>-z")
  out2
}

real_df = real_data(real_par)

plot0 <- function(df){
  ggplot(df )+geom_line(mapping = aes(x= times, y=x),color ='red')+
    geom_line(mapping = aes(x= times, y=y),color ='blue')+
    geom_line(mapping = aes(x= times, y=z),color ='green')+coord_fixed()
}

plot0(real_df)


sim_data <-  function(par,sd=0.01){
  r = c(par[1],par[4],par[7])
  A = matrix(c(0, par[2], par[3] ,par[5], 0, par[6], par[8], par[9], 0), nrow = 3, ncol = 3, byrow = T)
  x0 <- c(1, 1, 2)
  times = seq(4,7,length = nn)
  
  parameters <- list(r = r, A = A)
  # solve numerically
  out <- ode(y = x0, times = times, 
             func = GLV, parms = parameters)[,-1]
  out[2:nn,1] = out[2:nn,1] + rnorm(nn-1, sd = sd)
  out[2:nn,2] = out[2:nn,2] + rnorm(nn-1, sd = sd)
  out[2:nn,3] = out[2:nn,3] + rnorm(nn-1, sd = sd)
  out = data.frame(out)
  colnames(out) = c("x","y","z")
  out
}


df = replicate(30,sim_data(par),simplify = F)
times = seq(4,7,length = 30)

#df = df[[1]]; times = times


LV_3 <- function(df, times){
  power_equation <- function(x, power_par){ t(sapply(1:nrow(power_par),
                                                     function(c) power_par[c,1]*x^power_par[c,2] ) )}
  
  power_equation_fit <- function(data, X, n = 30){
    df_fit = t(sapply(1:nrow(data),function(c) power_equation_base(X, data[c,])))
    rownames(df_fit) = rownames(data)
    colnames(df_fit) = seq(min(X),max(X),length=n)
    return(df_fit)
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
  
  power_par = lapply(1:3, function(c) power_equation_base(y = df[,c], x = times))
  power_par = t(sapply(power_par, coef))
  
  #Time = times
  LotVmod <- function (Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      dx = x*Pars[1] + Pars[2]*y + Pars[3]*z
      dy = power_par[2,1]*power_par[2,2]*Time^(power_par[2,2]-1)
      dz = power_par[3,1]*power_par[3,2]*Time^(power_par[3,2]-1)
      return(list(c(dx, dy, dz)))
    })
  }
  #power_par[1,1]*times^(power_par[1,2]-1)
  
  f1 <- function(pars){
    #Pars <- rep(0.1,3)
    State <- c(x = df$x[1], y = df$y[1], z = df$z[1])
    Time <- times
    out <- as.data.frame(ode(func = LotVmod, y = State, parms = pars, times = Time ))[,-1]
    sse = crossprod(out[,1] - df[,1])
    sse
  }
  r1 = optim(rep(0.1,3), f1, method = "Nelder-Mead", control = list(maxit = 1e4, trace = T))
  #pars = r1$par
  
  LotVmod2 <- function (Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      dy = y*Pars[1] + Pars[2]*x + Pars[3]*z
      dx = power_par[1,1]*power_par[1,2]*Time^(power_par[1,2]-1)
      dz = power_par[3,1]*power_par[3,2]*Time^(power_par[3,2]-1)
      return(list(c(dx, dy, dz)))
    })
  }
  
  f2 <- function(pars){
    #Pars <- rep(0.1,3)
    State <- c(x = df$x[1], y = df$y[1], z = df$z[1])
    Time <- times
    out <- as.data.frame(ode(func = LotVmod2, y = State, parms = pars, times = Time))[,-1]
    sse = crossprod(out[,2] - df[,2])
    sse
  }
  r2 = optim(rep(0.1,3), f2, method = "Nelder-Mead", control = list(maxit = 1e4, trace = T))
  #pars = r1$par
  r2
  
  
  LotVmod3 <- function (Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      dz = z*Pars[1] + Pars[2]*x + Pars[3]*y
      dx = power_par[1,1]*power_par[1,2]*Time^(power_par[1,2]-1)
      dy = power_par[2,1]*power_par[2,2]*Time^(power_par[2,2]-1)
      return(list(c(dx, dy, dz)))
    })
  }
  
  f3 <- function(pars){
    #Pars <- rep(0.1,3)
    State <- c(x = df$x[1], y = df$y[1], z = df$z[1])
    Time <- times
    out <- as.data.frame(ode(func = LotVmod3, y = State, parms = pars, times = Time))[,-1]
    sse = crossprod(out[,3] - df[,3])
    sse
  }
  r3 = optim(rep(0.1,3), f3, method = "Nelder-Mead", control = list(maxit = 1e4, trace = T))
  #pars = r1$par
  r3
  par_est = list(par_df = data.frame(par1 = r1$par, par2 = r2$par, par3 = r3$par),
                 par = c(r1$par, r2$par, r3$par))
  return(par_est)
}
#par_est = LV_1(df = df[[1]], times = times)


core.number <- 16
cl <- makeCluster(getOption("cl.cores", core.number))
clusterEvalQ(cl, {library(deSolve)})
clusterExport(cl,c(c("GLV3","GLV","LV_3","df","times"),ls()),
              envir=environment())
par_est1 <- parLapply(cl=cl, 1:length(df), function(c)  LV_3(df = df[[c]], times))
stopCluster(cl)

#par_est1 = lapply(1:length(df), function(c) LV_1(df = df[[c]], times))
#par_est1= par_est2

LV_2 <- function(data,times,par0){
  GLV <- function(t, x, parameters){
    with(as.list(c(x, parameters)), {
      #x[x < 10^-8] <- 0 # prevent numerical problems
      dxdt <- x * r + A %*% x
      list(dxdt)
    })
  }
  
  ff <-  function(par, data, times){
    r = c(par[1],par[4],par[7])
    A = matrix(c(0, par[2], par[3] ,par[5], 0, par[6], par[8], par[9], 0), nrow = 3, ncol = 3, byrow = T)
    x0 <- as.numeric(data[1,])
    
    
    parameters <- list(r = r, A = A)
    # solve numerically
    out <- ode(y = x0, times = times, 
               func = GLV, parms = parameters)[,-1]
    
    sse = sum(sapply(1:3,function(c) crossprod(out[,c] - data[,c])))
    sse
  }
  
  s2 = optim(par0, ff, data = data, times = times,
             method = "L-BFGS-B", 
             lower = rep(-0.6,9),
             upper = rep(0.6,9),
             control = list(maxit = 1e5, trace = T))
  par = s2$par
  par_est = list(par_df = data.frame(par1 = par[1:3], par2 = par[4:6], par3 = par[7:9]),
                 par = par)
  return(par_est)
}


#par_est2 = lapply(1:length(df), function(c)LV_2(data = df[[c]], times = times, par0 = par_est1[[c]][[2]]))

core.number <- 16
cl <- makeCluster(getOption("cl.cores", core.number))
clusterEvalQ(cl, {library(deSolve)})
clusterExport(cl,c(c("GLV3","GLV","LV_2","df","times","par_est1"),ls()),
              envir=environment())
par_est2 <- parLapply(cl=cl, 1:length(df), function(c) LV_2(data = df[[c]], times = times, par0 = par_est1[[c]][[2]]))
stopCluster(cl)


par1 = colMeans(Reduce(rbind, lapply(par_est1, "[[", "par")))
par1sd = apply(Reduce(rbind, lapply(par_est1, "[[", "par")),2,sd)

par2 = colMeans(Reduce(rbind, lapply(par_est2, "[[", "par")))
par2sd = apply(Reduce(rbind, lapply(par_est2, "[[", "par")),2,sd)

real_par = c(real_par1,real_par2,real_par3)


par = colMeans(Reduce(rbind, lapply(par_est2, "[[", "par")))

darken <- function(color, factor=1.2){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

plot1 <- function(par, real_df){
  r = c(par[1],par[4],par[7])
  A = matrix(c(0, par[2], par[3] ,par[5], 0, par[6], par[8], par[9], 0), nrow = 3, ncol = 3, byrow = T)
  x0 <- c(1, 1, 2)
  times = seq(4,7,length = 30)
  
  parameters <- list(r = r, A = A)
  dep_list = list(c(2,3), c(1,3), c(1,2))
  n = 3
  x00 = c(x0,x0,rep(0,6))
  out2 = data.frame(ode(y = x00, times = times, 
                        func = GLV3, parms = parameters,
                        dep_list = dep_list))
  colnames(out2) = c("times","x","y","z","x_ind","y_ind","z_ind", 
                     "y>-x","z>-x",
                     "x>-y","z>-y",
                     "x>-z","y>-z")
  
  
  r1 = real_df[,c(1,2,5,8,9)]
  d1 = out2[,c(1,2,5,8,9)]
  rname1 = subset(d1, d1$times == 13)
  
  p1 = ggplot() + geom_line(r1, mapping = aes(x = times, y = x), color = darken("blue")) +
    geom_line(r1, mapping = aes(x = times, y = x_ind), color = darken("red")) +
    geom_line(r1, mapping = aes(x = times, y = r1$`y>-x`), color = darken("green")) +
    geom_line(r1, mapping = aes(x = times, y = r1$`z>-x`), color = darken("green")) +
    
    geom_line(d1, mapping = aes(x = times, y = x), color = darken("blue"),linetype=2) +
    geom_line(d1, mapping = aes(x = times, y = x_ind), color = darken("red"),linetype=2) +
    geom_line(d1, mapping = aes(x = times, y = d1$`y>-x`), color = darken("green"),linetype=2) +
    geom_line(d1, mapping = aes(x = times, y = d1$`z>-x`), color = darken("green"),linetype=2) +
    
    geom_text(rname1, mapping = aes(label = "y", x = rname1$times*1.01, y = rname1$`y>-x`), color = darken("green"))+
    geom_text(rname1, mapping = aes(label = "z", x = rname1$times*1.01, y = rname1$`z>-x`), color = darken("green"))+
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') + ylab("Gene Expression") + 
    ggtitle("x") + theme_bw() + scale_y_continuous(limits = c(-1,2)) + xlab(NULL) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0,0,0,0),"lines")) 
  
  
  r2 = real_df[,c(1,3,6,10,11)]
  d2 = out2[,c(1,3,6,10,11)]
  rname2 = subset(d2, d2$times == 13)
  
  p2 = ggplot() + geom_line(r2, mapping = aes(x = times, y = y), color = darken("blue")) +
    geom_line(r2, mapping = aes(x = times, y = y_ind), color = darken("red")) +
    geom_line(r2, mapping = aes(x = times, y = r2$`x>-y`), color = darken("green")) +
    geom_line(r2, mapping = aes(x = times, y = r2$`z>-y`), color = darken("green")) +
    
    geom_line(d2, mapping = aes(x = times, y = y), color = darken("blue"),linetype=2) +
    geom_line(d2, mapping = aes(x = times, y = y_ind), color = darken("red"),linetype=2) +
    geom_line(d2, mapping = aes(x = times, y = d2$`x>-y`), color = darken("green"),linetype=2) +
    geom_line(d2, mapping = aes(x = times, y = d2$`z>-y`), color = darken("green"),linetype=2) +
    
    geom_text(rname2, mapping = aes(label = "x", x = rname2$times*1.01, y = rname2$`x>-y`), color = darken("green"))+
    geom_text(rname2, mapping = aes(label = "z", x = rname2$times*1.01, y = rname2$`z>-y`), color = darken("green"))+
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') + ylab(NULL)  + xlab("Niche Index") +
    ggtitle("y")  + theme_bw() + scale_y_continuous(limits = c(-1,2))+ 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank(), 
          axis.ticks.length.y = unit(-0.1,"cm"),
          plot.margin = unit(c(0,0,0,0),"lines"))
  
  
  r3 = real_df[,c(1,4,7,12,13)]
  d3 = out2[,c(1,4,7,12,13)]
  rname3 = subset(d3, d3$times == 13)
  
  p3 = ggplot() + geom_line(r3, mapping = aes(x = times, y = z), color = darken("blue")) +
    geom_line(r3, mapping = aes(x = times, y = z_ind), color = darken("red")) +
    geom_line(r3, mapping = aes(x = times, y = r3$`x>-z`), color = darken("green")) +
    geom_line(r3, mapping = aes(x = times, y = r3$`y>-z`), color = darken("green")) +
    
    geom_line(d3, mapping = aes(x = times, y = z), color = darken("blue"),linetype=2) +
    geom_line(d3, mapping = aes(x = times, y = z_ind), color = darken("red"),linetype=2) +
    geom_line(d3, mapping = aes(x = times, y = d3$`x>-z`), color = darken("green"),linetype=2) +
    geom_line(d3, mapping = aes(x = times, y = d3$`y>-z`), color = darken("green"),linetype=2) +
    
    geom_text(rname3, mapping = aes(label = "x", x = rname3$times*1.01, y = rname3$`x>-z`), color = darken("green"))+
    geom_text(rname3, mapping = aes(label = "y", x = rname3$times*1.01, y = rname3$`y>-z`), color = darken("green"))+
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') + ylab(NULL)  + xlab(NULL) +
    ggtitle("z")  + theme_bw() + scale_y_continuous(limits = c(-1,2))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank(), 
          axis.ticks.length.y = unit(-0.1,"cm"),
          plot.margin = unit(c(0,0,0,0),"lines"))
  
  pp = p1+p2+p3
  pp
}

#p1 = plot1(colMeans(Reduce(rbind, lapply(par_est1, "[[", "par"))), real_df)
pA = plot1(colMeans(Reduce(rbind, lapply(par_est2, "[[", "par"))), real_df)



test_H0 <- function(par){
  mod0 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dx <-  a * x 
      dy <-  b * y
      dz <-  c * y
      list(c(dx,dy,dz))
    })
  }
  
  f0 <- function(init_pars){
    parameters <- c(a = init_pars[1], b = init_pars[2], c = init_pars[3])
    state <- c(x = real_df[1,2], y = real_df[1,3], z = real_df[1,4])
    times <- times
    out <- ode(y = state, times = times, func = mod0, parms = parameters)
    
    ssr = crossprod(real_df[,2]-out[,2])+
      crossprod(real_df[,3]-out[,3])+
      crossprod(real_df[,4]-out[,4])
    ssr
  }
  
  pars = optim(rep(0.1,3),f0)
  
  f0_fit <- function(pars){
    parameters <- c(a = pars[1], b = pars[2], c = pars[3])
    state <- c(x = real_df[1,2], y = real_df[1,3], z = real_df[1,4])
    times <- times
    out <- ode(y = state, times = times, func = mod0, parms = parameters)
    return(data.frame(out))
  }
  
  obs_fit = f0_fit(pars$par)
  
  
  r1 = real_df[,c(1,2,5,8,9)]
  
  p1 = ggplot() + geom_line(r1, mapping = aes(x = times, y = x), color = darken("blue")) +
    geom_line(r1, mapping = aes(x = times, y = x_ind), color = darken("red")) +
    geom_line(r1, mapping = aes(x = times, y = r1$`y>-x`), color = darken("green")) +
    geom_line(r1, mapping = aes(x = times, y = r1$`z>-x`), color = darken("green")) +
    
    geom_line(obs_fit, mapping = aes(time, x), color = "blue",linetype=2)+
    
    #geom_text(rname1, mapping = aes(label = "y", x = rname1$times*1.01, y = rname1$`y>-x`), color = darken("green"))+
    #geom_text(rname1, mapping = aes(label = "z", x = rname1$times*1.01, y = rname1$`z>-x`), color = darken("green"))+
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') + ylab("Gene Expression") + 
    ggtitle("x") + theme_bw() + scale_y_continuous(limits = c(-1,2)) + xlab(NULL) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0,0,0,0),"lines")) 
  
  
  r2 = real_df[,c(1,3,6,10,11)]
  
  p2 = ggplot() + geom_line(r2, mapping = aes(x = times, y = y), color = darken("blue")) +
    geom_line(r2, mapping = aes(x = times, y = y_ind), color = darken("red")) +
    geom_line(r2, mapping = aes(x = times, y = r2$`x>-y`), color = darken("green")) +
    geom_line(r2, mapping = aes(x = times, y = r2$`z>-y`), color = darken("green")) +
    
    geom_line(obs_fit, mapping = aes(time, y), color = "blue",linetype=2)+
    
    #geom_text(rname2, mapping = aes(label = "x", x = rname2$times*1.01, y = rname2$`x>-y`), color = darken("green"))+
    #geom_text(rname2, mapping = aes(label = "z", x = rname2$times*1.01, y = rname2$`z>-y`), color = darken("green"))+
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') + ylab(NULL)  + xlab("Niche Index") +
    ggtitle("y")  + theme_bw() + scale_y_continuous(limits = c(-1,2))+ 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank(), 
          axis.ticks.length.y = unit(-0.1,"cm"),
          plot.margin = unit(c(0,0,0,0),"lines"))
  
  
  r3 = real_df[,c(1,4,7,12,13)]
  
  p3 = ggplot() + geom_line(r3, mapping = aes(x = times, y = z), color = darken("blue")) +
    geom_line(r3, mapping = aes(x = times, y = z_ind), color = darken("red")) +
    geom_line(r3, mapping = aes(x = times, y = r3$`x>-z`), color = darken("green")) +
    geom_line(r3, mapping = aes(x = times, y = r3$`y>-z`), color = darken("green")) +
    
    geom_line(obs_fit, mapping = aes(time, z), color = "blue",linetype=2)+
    
    #geom_text(rname3, mapping = aes(label = "x", x = rname3$times*1.01, y = rname3$`x>-z`), color = darken("green"))+
    #geom_text(rname3, mapping = aes(label = "y", x = rname3$times*1.01, y = rname3$`y>-z`), color = darken("green"))+
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') + ylab(NULL)  + xlab(NULL) +
    ggtitle("z")  + theme_bw() + scale_y_continuous(limits = c(-1,2))+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank(), 
          axis.ticks.length.y = unit(-0.1,"cm"),
          plot.margin = unit(c(0,0,0,0),"lines"))
  
  pp = p1+p2+p3
  pp
  
}

pC = test_H0(real_par)


#2. data generate by nongame model------------------------------------
power_equation <- function(x, power_par){ t(sapply(1:nrow(power_par),
                                                   function(c) power_par[c,1]*x^power_par[c,2] ) )}

power_equation_fit <- function(data, X, n = 30){
  df_fit = t(sapply(1:nrow(data),function(c) power_equation_base(X, data[c,])))
  rownames(df_fit) = rownames(data)
  colnames(df_fit) = seq(min(X),max(X),length=n)
  return(df_fit)
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

power_par = lapply(1:3, function(c) power_equation_base(y = real_df[,c+1], x = times))
power_par = t(sapply(power_par, coef))

real_df2 <- data.frame(t(power_equation(times, power_par)))
colnames(real_df2) = c("x","y","z")

sim_data2 <-  function(sd=0.01){
  real_df2[2:nn,1] = real_df2[2:nn,1]+rnorm(nn-1, sd = sd)
  real_df2[2:nn,2] = real_df2[2:nn,2]+rnorm(nn-1, sd = sd)
  real_df2[2:nn,3] = real_df2[2:nn,3]+rnorm(nn-1, sd = sd)
  real_df2[1,1:3] = c(1, 1, 2) #same initial value
  real_df2
}


df2 = replicate(30,sim_data2(),simplify = F)



#analysis by Module 1

LV_11 <- function(data,times,par0){
  GLV <- function(t, x, parameters){
    with(as.list(c(x, parameters)), {
      #x[x < 10^-8] <- 0 # prevent numerical problems
      dxdt <- x * r + A %*% x
      list(dxdt)
    })
  }
  
  ff <-  function(par, data, times){
    r = c(par[1],par[4],par[7])
    A = matrix(c(0, par[2], par[3] ,par[5], 0, par[6], par[8], par[9], 0), nrow = 3, ncol = 3, byrow = T)
    x0 <- as.numeric(data[1,])
    
    
    parameters <- list(r = r, A = A)
    # solve numerically
    out <- ode(y = x0, times = times, 
               func = GLV, parms = parameters)[,-1]
    
    sse = sum(sapply(1:3,function(c) crossprod(out[,c] - data[,c])))
    sse
  }
  
  s2 = optim(par0, ff, data = data, times = times,
             method = "L-BFGS-B", 
             lower = rep(-0.6,9),
             upper = rep(0.6,9),
             control = list(maxit = 1e5, trace = T))
  par = s2$par
  par_est = list(par_df = data.frame(par1 = par[1:3], par2 = par[4:6], par3 = par[7:9]),
                 par = par)
  return(par_est)
}

core.number <- 16
cl <- makeCluster(getOption("cl.cores", core.number))
clusterEvalQ(cl, {library(deSolve)})
clusterExport(cl,c(c("GLV3","GLV","LV_11","df2","times","par_est1"),ls()),
              envir=environment())
par_est3 <- parLapply(cl=cl, 1:length(df2), function(c) LV_11(data = df2[[c]], times = times, par0 = par_est1[[c]][[2]]))
stopCluster(cl)


par3 = colMeans(Reduce(rbind, lapply(par_est3, "[[", "par")))
par3sd = apply(Reduce(rbind, lapply(par_est3, "[[", "par")),2,sd)


plot2 <- function(par, real_df){
  r = c(par[1],par[4],par[7])
  A = matrix(c(0, par[2], par[3] ,par[5], 0, par[6], par[8], par[9], 0), nrow = 3, ncol = 3, byrow = T)
  x0 <- c(1, 1, 2)
  times = seq(4,7,length = 30)
  
  parameters <- list(r = r, A = A)
  dep_list = list(c(2,3), c(1,3), c(1,2))
  n = 3
  x00 = c(x0,x0,rep(0,6))
  out2 = data.frame(ode(y = x00, times = times, 
                        func = GLV3, parms = parameters,
                        dep_list = dep_list))
  colnames(out2) = c("times","x","y","z","x_ind","y_ind","z_ind", 
                     "y>-x","z>-x",
                     "x>-y","z>-y",
                     "x>-z","y>-z")
  
  real_df = cbind(times,real_df)
  colnames(real_df)[1]="times"
  r1 = real_df[,c(1,2)]
  d1 = out2[,c(1,2,5,8,9)]
  rname1 = subset(d1, d1$times == 13)
  
  p1 = ggplot() + geom_line(r1, mapping = aes(x = times, y = x), color = darken("blue")) +
    #geom_line(r1, mapping = aes(x = times, y = x_ind), color = darken("red")) #+
    #geom_line(r1, mapping = aes(x = times, y = r1$`y>-x`), color = darken("green")) +
    #geom_line(r1, mapping = aes(x = times, y = r1$`z>-x`), color = darken("green")) +
    
    geom_line(d1, mapping = aes(x = times, y = x), color = darken("blue"),linetype=2) +
    geom_line(d1, mapping = aes(x = times, y = x_ind), color = darken("red"),linetype=2) +
    geom_line(d1, mapping = aes(x = times, y = d1$`y>-x`), color = darken("green"),linetype=2) +
    geom_line(d1, mapping = aes(x = times, y = d1$`z>-x`), color = darken("green"),linetype=2) +
    
    geom_text(rname1, mapping = aes(label = "y", x = rname1$times*1.01, y = rname1$`y>-x`), color = darken("green"))+
    geom_text(rname1, mapping = aes(label = "z", x = rname1$times*1.01, y = rname1$`z>-x`), color = darken("green"))+
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') + ylab("Gene Expression") + 
    ggtitle("x") + theme_bw() + scale_y_continuous(limits = c(-1,2)) + xlab(NULL) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0,0,0,0),"lines")) 
  
  
  r2 =  real_df[,c(1,3)]
  d2 = out2[,c(1,3,6,10,11)]
  rname2 = subset(d2, d2$times == 13)
  
  p2 = ggplot() + geom_line(r2, mapping = aes(x = times, y = y), color = darken("blue")) +
    #geom_line(r2, mapping = aes(x = times, y = y_ind), color = darken("red")) +
    #geom_line(r2, mapping = aes(x = times, y = r2$`x>-y`), color = darken("green")) +
    #geom_line(r2, mapping = aes(x = times, y = r2$`z>-y`), color = darken("green")) +
    
    geom_line(d2, mapping = aes(x = times, y = y), color = darken("blue"),linetype=2) +
    geom_line(d2, mapping = aes(x = times, y = y_ind), color = darken("red"),linetype=2) +
    geom_line(d2, mapping = aes(x = times, y = d2$`x>-y`), color = darken("green"),linetype=2) +
    geom_line(d2, mapping = aes(x = times, y = d2$`z>-y`), color = darken("green"),linetype=2) +
    
    geom_text(rname2, mapping = aes(label = "x", x = rname2$times*1.01, y = rname2$`x>-y`), color = darken("green"))+
    geom_text(rname2, mapping = aes(label = "z", x = rname2$times*1.01, y = rname2$`z>-y`), color = darken("green"))+
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') + ylab(NULL)  + xlab("Niche Index") +
    ggtitle("y")  + theme_bw() + scale_y_continuous(limits = c(-1,2)) + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank(), 
          axis.ticks.length.y = unit(-0.1,"cm"),
          plot.margin = unit(c(0,0,0,0),"lines"))
  
  
  r3 = real_df[,c(1,4)]
  d3 = out2[,c(1,4,7,12,13)]
  rname3 = subset(d3, d3$times == 13)
  
  p3 = ggplot() + geom_line(r3, mapping = aes(x = times, y = z), color = darken("blue")) +
    #geom_line(r3, mapping = aes(x = times, y = z_ind), color = darken("red")) +
    #geom_line(r3, mapping = aes(x = times, y = r3$`x>-z`), color = darken("green")) +
    #geom_line(r3, mapping = aes(x = times, y = r3$`y>-z`), color = darken("green")) +
    
    geom_line(d3, mapping = aes(x = times, y = z), color = darken("blue"),linetype=2) +
    geom_line(d3, mapping = aes(x = times, y = z_ind), color = darken("red"),linetype=2) +
    geom_line(d3, mapping = aes(x = times, y = d3$`x>-z`), color = darken("green"),linetype=2) +
    geom_line(d3, mapping = aes(x = times, y = d3$`y>-z`), color = darken("green"),linetype=2) +
    
    geom_text(rname3, mapping = aes(label = "x", x = rname3$times*1.01, y = rname3$`x>-z`), color = darken("green"))+
    geom_text(rname3, mapping = aes(label = "y", x = rname3$times*1.01, y = rname3$`y>-z`), color = darken("green"))+
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') + ylab(NULL)  + xlab(NULL) +
    ggtitle("z")  + theme_bw()+ scale_y_continuous(limits = c(-1,2))  +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank(), 
          axis.ticks.length.y = unit(-0.1,"cm"),
          plot.margin = unit(c(0,0,0,0),"lines"))
  
  pp = p1+p2+p3
  pp
}

#p1 = plot1(colMeans(Reduce(rbind, lapply(par_est1, "[[", "par"))), real_df)
pB = plot2(colMeans(Reduce(rbind, lapply(par_est3, "[[", "par"))), real_df2)



plot3 <- function(par, real_df){
  real_df = cbind(times,real_df)
  colnames(real_df)[1]="times"
  mod0 <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dx <-  a * x 
      dy <-  b * y
      dz <-  c * y
      list(c(dx,dy,dz))
    })
  }
  
  f0 <- function(init_pars){
    parameters <- c(a = init_pars[1], b = init_pars[2], c = init_pars[3])
    state <- c(x = real_df[1,2], y = real_df[1,3], z = real_df[1,4])
    times <- times
    out <- ode(y = state, times = times, func = mod0, parms = parameters)
    
    ssr = crossprod(real_df[,2]-out[,2])+
      crossprod(real_df[,3]-out[,3])+
      crossprod(real_df[,4]-out[,4])
    ssr
  }
  
  pars = optim(rep(0.1,3),f0)
  
  f0_fit <- function(pars){
    parameters <- c(a = pars[1], b = pars[2], c = pars[3])
    state <- c(x = real_df[1,2], y = real_df[1,3], z = real_df[1,4])
    times <- times
    out <- ode(y = state, times = times, func = mod0, parms = parameters)
    return(data.frame(out))
  }
  
  obs_fit = f0_fit(pars$par)
  
  
  p1 = ggplot() + geom_line(real_df, mapping = aes(x = times, y = x), color = darken("blue")) +
    geom_line(obs_fit, mapping = aes(time, x), color = darken("blue"),linetype=2) +
    
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') + ylab("Gene Expression") + 
    ggtitle("x") + theme_bw() + scale_y_continuous(limits = c(-1,2)) + xlab(NULL) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0,0,0,0),"lines")) 
  
  
  p2 = ggplot() + geom_line(real_df, mapping = aes(x = times, y = y), color = darken("blue")) +
    #geom_line(r2, mapping = aes(x = times, y = y_ind), color = darken("red")) +
    #geom_line(r2, mapping = aes(x = times, y = r2$`x>-y`), color = darken("green")) +
    #geom_line(r2, mapping = aes(x = times, y = r2$`z>-y`), color = darken("green")) +
    geom_line(obs_fit, mapping = aes(time, y), color = darken("blue"),linetype=2) +
    
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') + ylab(NULL)  + xlab("Niche Index") +
    ggtitle("y")  + theme_bw() + scale_y_continuous(limits = c(-1,2)) + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank(), 
          axis.ticks.length.y = unit(-0.1,"cm"),
          plot.margin = unit(c(0,0,0,0),"lines"))
  
  
  p3 = ggplot() + geom_line(real_df, mapping = aes(x = times, y = z), color = darken("blue")) +
    geom_line(obs_fit, mapping = aes(time, z), color = darken("blue"),linetype=2) +
    
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') + ylab(NULL)  + xlab(NULL) +
    ggtitle("z")  + theme_bw()+ scale_y_continuous(limits = c(-1,2))  +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank(), 
          axis.ticks.length.y = unit(-0.1,"cm"),
          plot.margin = unit(c(0,0,0,0),"lines"))
  
  pp = p1+p2+p3
  pp
}

pD = plot3(real_par,real_df2)


library(cowplot)
pp=plot_grid(pA,pB,pC,pD, labels = c('A', 'B','C','D'))
ggsave("FigS6.1.pdf",pp,width = 20, height = 8)
