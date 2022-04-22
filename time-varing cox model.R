library("MASS")
library(tidyverse)
library("extraDistr")
library(nleqslv)
library(nloptr)
library(splines)
library(gaussquad)


HPP_gen <- function( censor, lambda.bar) {
  
    nn <- rpois(1,  lambda.bar)+1
     tt <- sort(runif(nn, min = 0, max = 1))
    tt <- tt[which(tt<=censor)]
  tt
}

kerfun <- function(xx){
  pmax((1-xx^2)*0.75,0)
}

lqrule64 <- legendre.quadrature.rules(64)[[64]]

##################generate data############################
meanfunction<-function(t){-1-2*(t-0.5)^2}
simAsytransdata <-
  function(
           alpha,
           beta,
           cen,
           nstep  ) {
    #browser()
    # 1. Generate censoring times
    cen=runif(1,cen,1.5)
    cen <- min(cen,1)  
    
    # 2. Generate Z(t) and h(t) as a step function
    
    Sigmat_z   <-   exp(-abs(outer(1:nstep, 1:nstep, "-")) / nstep)
   
    z   <-  (c(mvrnorm(
      1, meanfunction(c(1:20)/20), Sigmat_z 
    )))
    
     
     
    left_time_points <- (0:(nstep - 1)) / nstep
    
    # 3. Generate observation time
   
    z_fun <- stepfun(left_time_points, c(0,z  ))
    h_fun <- function(x) {
      beta(x) * z_fun(x)
    }
    
    # 4. Generate failure time
    
    lam_fun <- function(tt)
      alpha(tt) * exp(h_fun(tt))
    
    u <- runif(1)
    
    fail_time <- nleqslv(0, function(ttt)
      legendre.quadrature(lam_fun,
                          lower = 0,
                          upper = ttt,
                          lqrule64) + log(u))$x
    
    
    
    X <- min(fail_time, cen)
    
    # 1. Generate R_ik
    
    obs_times <- HPP_gen(
                 cen, 5)
    
    if (length(obs_times) == 0)
      obs_times <- cen
    
    covariates_obscov <-z_fun(obs_times)
    
    # Results
    return(
      tibble(
        X = X,
        delta = fail_time < cen,
        covariates = covariates_obscov,
        obs_times = obs_times,
        censoring = cen
      ) # %>% filter(X >= obs_times)
    )
    
  }


##################estimating beta_0(t)############################
estproc_Cox <- function(data, n,s, h1,h2, pl = 0) {
  X=data$X
  id <- data$id
  covariates <- data$covariates 
  obs_times <- data$obs_times
  delta <- data$delta
  kerval <- kerfun((X - s) / h1) * kerfun((obs_times-s)/h2) / h1 /h2 

  
  # Estimating equation and estimation
  
  Delta<-outer(X,X,"<=")
  kerval_XX<- kerfun((X-s)/h1) %*%t(kerfun((obs_times-s)/h2)  )/h1/h2
  
  estequ <- function(beta) {
    
    expbeta <- as.vector(exp(beta*covariates  ))
    S0_XX <-  t(kerval_XX * Delta)* expbeta    
    Zbar_XX <-  t(S0_XX) %*% covariates/  colSums(S0_XX)
    Zbar_XX[is.na(Zbar_XX)] <- 0 
    res <- sum(delta*kerval*(covariates-Zbar_XX))
    res
  }
  
  estres <- nleqslv(0,fn=estequ)
  
  beta_est <- estres$x
  expbeta <- as.vector(exp(beta_est*covariates))
  S0_XX <-  t(kerval_XX * Delta)* expbeta 
  Zbar_XX <-   t(S0_XX) %*% covariates/  colSums(S0_XX)
  Zbar_XX[is.na(Zbar_XX)] <- 0
 # dZbar_XX <-   (t(S0_XX) %*% covariates^2*colSums(S0_XX)-(t(S0_XX) %*% covariates)^2)/  (colSums(S0_XX))^2
  dZbar_XX=  ( t(S0_XX) %*% (covariates^2)*colSums(S0_XX)/colSums(S0_XX)^2-( t(S0_XX) %*% covariates)^2/colSums(S0_XX)^2)
  dZbar_XX[is.na(dZbar_XX)] <- 0 
  
  B_indi <- tibble(B_indi_inner = delta*kerval*(covariates-Zbar_XX),subid= as.integer(id)) %>% group_by(subid) %>% summarise(B=t(sum(B_indi_inner)))
  B <- sum(B_indi$B^2)/n^2
  A_indi <- tibble(A_indi_inner = -delta*kerval*dZbar_XX,subid= as.integer(id)) %>% group_by(subid) %>% summarise(A=t(sum(A_indi_inner)))
  A <- sum(A_indi$A)/n

 # se=sqrt(sum((kerval*(covariates-Zbar_XX))^2)/sum(kerval*dZbar_XX)^2) 
  
  sigma=ifelse(A*B==0,0,sqrt(solve(A) %*% B %*% solve(A))) 
  #list(est=beta_est,se=sigma)
  M=5000
  xi=matrix(rexp(M*n,1)-1,n)
  #xi=xi[as.numeric(id),]
  res=colSums(c(B_indi$B)*xi) 
  
  tilde_s=abs(c(1/sigma/A )*(res))/n

  
  list(est=beta_est,sigma=sigma,tilde_s=tilde_s)
  
}
estproc_Cox_noscb  <- function(data, n,s, h1,h2, pl = 0) {
  X=data$X
  id <- data$id
  covariates <- data$covariates 
  obs_times <- data$obs_times
  delta <- data$delta
  kerval <- kerfun((X - s) / h1) * kerfun((obs_times-s)/h2) / h1 /h2 
  
  
  # Estimating equation and estimation
  
  Delta<-outer(X,X,"<=")
  kerval_XX<- kerfun((X-s)/h1) %*%t(kerfun((obs_times-s)/h2)  )/h1/h2
  
  estequ <- function(beta) {
    
    expbeta <- as.vector(exp(beta*covariates  ))
    S0_XX <-  t(kerval_XX * Delta)* expbeta    
    Zbar_XX <-  t(S0_XX) %*% covariates/  colSums(S0_XX)
    Zbar_XX[is.na(Zbar_XX)] <- 0 
    res <- sum(delta*kerval*(covariates-Zbar_XX))
    res
  }
  
  estres <- nleqslv(0,fn=estequ)
  
  beta_est <- estres$x
  expbeta <- as.vector(exp(beta_est*covariates))
  S0_XX <-  t(kerval_XX * Delta)* expbeta 
  Zbar_XX <-   t(S0_XX) %*% covariates/  colSums(S0_XX)
  Zbar_XX[is.na(Zbar_XX)] <- 0
  # dZbar_XX <-   (t(S0_XX) %*% covariates^2*colSums(S0_XX)-(t(S0_XX) %*% covariates)^2)/  (colSums(S0_XX))^2
  dZbar_XX=  ( t(S0_XX) %*% (covariates^2)*colSums(S0_XX)/colSums(S0_XX)^2-( t(S0_XX) %*% covariates)^2/colSums(S0_XX)^2)
  dZbar_XX[is.na(dZbar_XX)] <- 0 
  
  B_indi <- tibble(B_indi_inner = delta*kerval*(covariates-Zbar_XX),subid= as.integer(id)) %>% group_by(subid) %>% summarise(B=t(sum(B_indi_inner)))
  B <- sum(B_indi$B^2)/n^2
  A_indi <- tibble(A_indi_inner = -delta*kerval*dZbar_XX,subid= as.integer(id)) %>% group_by(subid) %>% summarise(A=t(sum(A_indi_inner)))
  A <- sum(A_indi$A)/n
  
  sigma=ifelse(A*B==0,0,sqrt(solve(A) %*% B %*% solve(A))) 

  list(est=beta_est,sigma=sigma)
  
}
estproc_Cox_test <- function(data, n,s, h1,h2, pl = 0) {
  X=data$X
  id <- data$id
  covariates <- data$covariates 
  obs_times <- data$obs_times
  delta <- data$delta
  kerval <- kerfun((X - s) / h1) * kerfun((obs_times-s)/h2) / h1 /h2 
  
  
  # Estimating equation and estimation
  
  Delta<-outer(X,X,"<=")
  kerval_XX<- kerfun((X-s)/h1) %*%t(kerfun((obs_times-s)/h2)  )/h1/h2
  
  estequ <- function(beta) {
    
    expbeta <- as.vector(exp(beta*covariates  ))
    S0_XX <-  t(kerval_XX * Delta)* expbeta    
    Zbar_XX <-  t(S0_XX) %*% covariates/  colSums(S0_XX)
    Zbar_XX[is.na(Zbar_XX)] <- 0 
    res <- sum(delta*kerval*(covariates-Zbar_XX))
    res
  }
  
  estres <- nleqslv(0,fn=estequ)
  
  beta_est <- estres$x
  
  
  list(est=beta_est)
  
}

##################### bandwidth selection ##################
hcv<-function(simdata,nn){
  test_idk <- lapply(split(sample(1:nn,nn),rep(1:2,nn/2)),sort)
  testnn=10
  hmin<-nn^(-0.45)
  hmax<-nn^(-0.35)
  hn <- exp(seq(log(hmin),log(hmax),length.out=testnn+1)[-1])
  ss<-c(1:10)/10*(1-2*hmax)+hmax
  bias1=hmax*hn
  bias2=hn^2
  res <-foreach(s=ss)%do%{
    reshh=foreach(hh=hn) %do% {
      beta1=estproc_Cox_test(simdata %>% filter(!(as.numeric(id) %in% test_idk$`1`)) ,
                                    nn/2,s,hmax,hh)$est
       
      beta2=estproc_Cox_test(simdata %>% filter(!(as.numeric(id) %in% test_idk$`2`)),
                                    nn/2,s,hmax,hh)$est
      sbeta=estproc_Cox_test(simdata,
                                    nn,s ,hmax,hh)$est
      list(beta1,beta2,sbeta)
    }
    reshh=reshh %>%
      lapply(function(xx) tibble(beta1=xx[[1]],beta2=xx[[2]],sbeta=xx[[3]])) %>% 
      bind_rows(.id="rep") %>% 
      mutate(var=(beta1-beta2)^2/4)
    C=lm(reshh$sbeta~bias1 +bias2 )$coefficients[2:3]
    return( (C%*%t(matrix(c(bias1,bias2),10)))^2 +reshh$var)
  }
  
  res1=data.frame(t(sapply(res,c)))
  res=hn[which.min(colSums(res1))]
  return(res)  
  
}


 
 


library(tidyverse)
library(caret)
library(foreach)
library(doParallel)
cores=detectCores()
cl <- makeCluster(cores)
doParallel::registerDoParallel(cl)
 
nn=400
s=.2
beta=function(t)  {0.5*sin(2*pi*t)   }
estres_all_h <- foreach(ii = 1:1000,.packages = c("tidyverse","extraDistr","MASS","nleqslv","nloptr","splines","gaussquad","foreach")) %dopar% {
             
           print(ii)
          
          simdata <-
            replicate(
              nn,
              simAsytransdata(
                alpha = function(tt) 2 +0.1*tt, 
                beta, 
                cen=0.7 #tuning the censoring rates 
                , 
                nstep = 20
              ),
              simplify = F
            ) %>% bind_rows(.id = "id")
           
           #  simdata %>% group_by(id) %>% filter(row_number() == 1)   %>% pull(delta) %>% mean()
           
           
     
            temp <- estproc_Cox(data=simdata,n=nn,s,h1=nn^(-0.35),h2=nn^(-0.4))
 

          list(temp)

}

truevalue=beta(s)  

estres_all_h %>%
  lapply(function(xx) tibble(est=xx[[1]]$est,se=xx[[1]]$sigma)) %>% 
  bind_rows(.id="rep") %>% 
  mutate(lb=est-qnorm(0.975)*se,
         ub=est+qnorm(0.975)*se,
         CP = (lb < truevalue)*(ub>truevalue)) %>% 
  summarise(bias=mean(est)-truevalue,sd=sd(est),se=mean(se),CP=mean(CP))

sd((estres_all_h  %>%
      lapply(function(xx) tibble(est=xx[[1]]$est,se=xx[[1]]$sigma)) %>% 
      bind_rows(.id="rep") %>% 
      mutate(lb=est-qnorm(0.975)*se,
             ub=est+qnorm(0.975)*se,
             CP = (lb < truevalue)*(ub>truevalue)))$est)
mean((estres_all_h  %>%
        lapply(function(xx) tibble(est=xx[[1]]$est,se=xx[[1]]$sigma)) %>% 
        bind_rows(.id="rep") %>% 
        mutate(lb=est-qnorm(0.975)*se,
               ub=est+qnorm(0.975)*se,
               CP = (lb < truevalue)*(ub>truevalue)))$se)

################################auto###########################################
 estres_all_h <- foreach(ii = 1:1000,.packages = c("tidyverse","extraDistr","MASS","nleqslv","nloptr","splines","gaussquad","foreach")) %dopar% {
  
  print(ii)
  
  simdata <-
    replicate(
      nn,
      simAsytransdata(
        alpha = function(tt) 2+0.1*tt, 
        beta, 
        cen=.7 
        , 
        nstep = 20
      ),
      simplify = F
    ) %>% bind_rows(.id = "id")
  
  hh=hcv(simdata,nn)
  
  
  temp <- estproc_Cox(data=simdata,n=nn,s,h1=nn^(-0.35),h2=hh)
  
  
  list(temp)
  
}

truevalue=beta(s)  

estres_all_h %>%
  lapply(function(xx) tibble(est=xx[[1]]$est,se=xx[[1]]$sigma)) %>% 
  bind_rows(.id="rep") %>% 
  mutate(lb=est-qnorm(0.975)*se,
         ub=est+qnorm(0.975)*se,
         CP = (lb < truevalue)*(ub>truevalue)) %>% 
  summarise(bias=mean(est)-truevalue,sd=sd(est),se=mean(se),CP=mean(CP))

sd((estres_all_h  %>%
      lapply(function(xx) tibble(est=xx[[1]]$est,se=xx[[1]]$sigma)) %>% 
      bind_rows(.id="rep") %>% 
      mutate(lb=est-qnorm(0.975)*se,
             ub=est+qnorm(0.975)*se,
             CP = (lb < truevalue)*(ub>truevalue)))$est)
mean((estres_all_h  %>%
        lapply(function(xx) tibble(est=xx[[1]]$est,se=xx[[1]]$sigma)) %>% 
        bind_rows(.id="rep") %>% 
        mutate(lb=est-qnorm(0.975)*se,
               ub=est+qnorm(0.975)*se,
               CP = (lb < truevalue)*(ub>truevalue)))$se)




############################### finding the threshold ###################################
nn=400
from=nn^(-0.35)
to=1-nn^(-0.35)
points=seq(from,to,by=(to-from)/49)
res=foreach(s=points) %do%{
  estres_all_h <- foreach(ii = 1:50 ,.packages = c("tidyverse","extraDistr","MASS","nleqslv","nloptr","splines","gaussquad","foreach")) %dopar% {
    
    print(ii)
    
    simdata <-
      replicate(
        nn,
        simAsytransdata(
          alpha = function(tt) 2+0.1*tt,
          beta,
          cen= .7,
          nstep = 20
        ),
        simplify = F
      ) %>% bind_rows(.id = "id")
    
    estres_ori <-  estproc_Cox(
      data = simdata,
      nn,
      s,
      h1=nn^(-0.35),
      h2=nn^(-0.45)
    )
    
    
    list(estres_ori)
    
  }
  res=estres_all_h %>%
    lapply(function(xx) tibble(est=xx[[1]]$est,sigma=xx[[1]]$sigma) )%>% 
    bind_rows(.id="rep")  %>% summarise(est=mean(est),sigma=mean(sigma))
  scb=estres_all_h %>%
    lapply(function(xx) tibble(est=xx[[1]]$tilde_s)) %>%
    bind_rows(.id="rep")  
  list(res,scb)
}



cc=res %>% lapply(function(xx) tibble(est=xx[[2]]$est)) %>% bind_cols(.id = NULL)  
M=5000
id=NULL
for(i in 1:1000){
  id=c(id,rep(i,M))
}
c_alpha=
  tibble(ss=apply(cc,1,max)) %>% mutate(id=id) %>% group_by(id)  %>%  summarise(q=quantile(ss,0.95)) 
meanq=mean(c_alpha$q)
 

##################################[ scb cp ]############################################################################

res=foreach(s=points) %do%{
  estres_all_h <- foreach(ii = 1:1000,.packages = c("tidyverse","extraDistr","MASS","nleqslv","nloptr","splines","gaussquad","foreach")) %dopar% {
    
    print(ii)
    
    simdata <-
      replicate(
        nn,
        simAsytransdata(
          alpha = function(tt) 2+0.1*tt,
          beta ,
          cen= 0.7,
          nstep = 20
        ),
        simplify = F
      ) %>% bind_rows(.id = "id")
    
    estres_ori <-  estproc_Cox_noscb(
      data = simdata,
      nn,
      s,
      h1=nn^(-0.35),
      h2=nn^(-0.45)
    )
    
    
    list(estres_ori)
    
  }
  res=estres_all_h %>%
    lapply(function(xx) tibble(est=xx[[1]]$est,sigma=xx[[1]]$sigma)) %>% bind_rows() 
  
  list(res)
}
testcoxn400bd3545cen20=res
save(testcoxn400bd3545cen20,file="testcoxn400bd3545cen20.RData")



true=NULL
for(i in 1:length(points))
{true=c(true,rep(0.5*sin(2*pi*points[i]),1000))}
est=testcoxn400bd3545cen20 %>%
  lapply(function(xx) tibble(est=xx[[1]]$est,sigma=xx[[1]]$sigma,lp=xx[[1]]$est-xx[[1]]$sigma*meanq,
                             up=xx[[1]]$est+xx[[1]]$sigma*meanq)) %>% 
  bind_rows()  %>% mutate(i=(up>true)*(lp<true)) %>% pull(i )   

sum(apply(matrix(est,1000),1,min))/1000

est=testcoxn400bd3545cen20 %>%
  lapply(function(xx) tibble(est=xx[[1]]$est,sigma=xx[[1]]$sigma,lp=xx[[1]]$est-xx[[1]]$sigma*1.96,
                             up=xx[[1]]$est+xx[[1]]$sigma*1.96)) %>% 
  bind_rows()  %>% mutate(i=(up>true)*(lp<true)) %>% pull(i )   

sum(apply(matrix(est,1000),1,min))/1000



     