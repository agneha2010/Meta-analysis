library(openxlsx)
library(metafor)
library(lme4)
library(reshape2)
library(dplyr)
library(locfit)
library(mvtnorm)
library(maxLik)
library(dplyr)

mi.data = read.xlsx("~meta-analysis/real_data/glm/data.xlsx")
k = nrow(mi.data)
mi.data$alive.control = mi.data$total.control-mi.data$death.control
mi.data$alive.treatment = mi.data$total.treatment-mi.data$death.treatment

data = list()
sigma.beta.hat.k.m <- c(0) ; beta.hat.k.m <- c(0); alpha.hat.k.m = c(0)
beta.hat.k <- list() ; sigma.beta.hat.k <-list()

for (i in 1:k){
  #i=1
  x = c(rep(1,mi.data$death.control[i]),rep(2,mi.data$death.treatment[i]),
        rep(1,mi.data$alive.control[i]),rep(2,mi.data$alive.treatment[i]))
  y = c(rep(1,mi.data$death.control[i]+mi.data$death.treatment[i]),
        rep(0,mi.data$alive.control[i]+mi.data$alive.treatment[i]))
  data[[i]] = data.frame(x=x,y=y)
  mylogit <- glm(y~x, data = data[[i]], family = "binomial")
  beta.hat.k[[i]] <- mylogit$coefficients
  sigma.beta.hat.k[[i]] <- vcov(mylogit)
  alpha.hat.k.m[i] <- beta.hat.k[[i]][[1]]
  beta.hat.k.m[i] <- beta.hat.k[[i]][[2]]
  sigma.beta.hat.k.m[i] <- sigma.beta.hat.k[[i]][2,2]
}
ini.sigma.alpha = var(alpha.hat.k.m)
ini.sigma.beta = var(beta.hat.k.m)
ini.beta = mean(beta.hat.k.m)
ini.alpha = mean(alpha.hat.k.m)

comb.est.per.ipd <- function(ipd.seq){
  ma.seq = which(!1:k %in% ipd.seq)
  k1 <- length(ipd.seq)
  k2 <- length(ma.seq)
  if(k1>0){
    beta.hat.k.ipd <- beta.hat.k[ipd.seq]
    sigma.beta.hat.k.ipd <- sigma.beta.hat.k[ipd.seq]
  }
  if(k2>0){
    beta.hat.k.meta <- beta.hat.k.m[ma.seq]
    sigma.beta.hat.k.meta <- sigma.beta.hat.k.m[ma.seq]
  }
 
  ####---- Using the closed form expressions
  theta.f = function(sigma.alpha,sigma.beta){
    s1 <- list() ; s2 <- list() ; m <- list() ; l <- list()
    if(k1>0){
      for (i in 1:k1){
        m[[i]] <- solve(matrix(c(sigma.alpha,0,0,sigma.beta),nrow=2) + as.matrix(sigma.beta.hat.k.ipd[[i]]))
        s1[[i]] <- m[[i]]%*%beta.hat.k.ipd[[i]]
      }
    }
    if(k2>0){
      for (i in 1:k2){
        l[[i]] <- matrix(c(0,0,0,(sigma.beta + sigma.beta.hat.k.meta[i])^(-1)),nrow=2)
        s2[[i]] <- l[[i]]%*%c(0,beta.hat.k.meta[i])
      }
    }
    if(k2==0){
      var = solve(Reduce("+",m))
      est = var%*%(Reduce("+",s1))
    } else if(k1==0){
      var = matrix(c(NA,NA,NA,solve(Reduce("+",l)[2,2])),nrow=2)
      est = c(NA,solve(Reduce("+",l)[2,2])%*%(Reduce("+",s2)[2]))
    }else{
      tmp1 = if (length(Reduce("+",m))==0) 0 else Reduce("+",m) + if (length(Reduce("+",l))==0) 0 else Reduce("+",l)
      tmp2 <- if (length(Reduce("+",s1))==0) 0 else Reduce("+",s1) + if (length(Reduce("+",s2))==0) 0 else Reduce("+",s2)
      var <- solve(tmp1)
      est <- var%*%tmp2
    }
    return(list(est,var))
  }
  
  logLikFun.1 <- function(par){
    #par = c(2,2)
    sigma.alpha <- par[1]
    sigma.beta <- par[2]
    theta = theta.f(sigma.alpha,sigma.beta)
    alpha = theta[[1]][1]
    beta = theta[[1]][2]
    s1 <- c(0) ; s2 <- c(0)
    if(k1>0){
      for (i in 1:k1){
        m <- matrix(c(sigma.alpha,0,0,sigma.beta),nrow=2) + as.matrix(sigma.beta.hat.k.ipd[[i]])
        s1[i] <- dmvnorm(beta.hat.k.ipd[[i]], mean = c(alpha,beta), sigma = m, log = TRUE)
      }
    }
    if(k2>0){
      for (i in 1:k2){
        l <- sqrt(sigma.beta + sigma.beta.hat.k.meta[i])
        s2[i] <- dnorm(beta.hat.k.meta[i], mean = beta, sd = l, log = TRUE)
      }
    }
    return(sum(s1,na.rm = T) + sum(s2, na.rm = T))
  }
  constraints <- list(ineqA=matrix(c(1,0,0,1),2,2,byrow=T), ineqB=c(0,0))
  if(k1>0){
    mle.1 <- maxLik(logLik=logLikFun.1,start=c(ini.sigma.alpha,ini.sigma.beta),constraints=constraints,method="BFGS")
  }else{
    mle.1 <- maxLik(logLik=logLikFun.1,start=c(sigma.alpha=ini.sigma.alpha,sigma.beta=ini.sigma.beta),
                    fixed=c("sigma.alpha"),constraints=constraints,method="BFGS")
  }
  gg.1 <- summary(mle.1)
  est.1 <- gg.1$estimate[,1]
  sd.1 <- gg.1$estimate[,2]
  if(k1==0){
    est.1[1] = NA
    sd.1[1] = NA
  }
  est.2 = theta.f(est.1[1],est.1[2])
  est.final = c(est.2[[1]][1],est.2[[1]][2],est.1[1],est.1[2])
  sd.final = c(sqrt(est.2[[2]][1,1]),sqrt(est.2[[2]][2,2]),sd.1[1],sd.1[2])
  mm = c(est.final,sd.final)
  return(mm)
}

prop.ipd <- c(0,0.2,0.4,0.6,0.8,1)
ipd.seq.list = list(NULL,c(1,2,17,19),c(1:3,17:22),1:13,1:18,1:22)
output <- round(as.data.frame(mapply(comb.est.per.ipd,ipd.seq.list)),4)
colnames(output) = paste0("prop.ipd_",prop.ipd)
rownames(output) = c("alpha","beta","sigma.alpha","sigma.beta",
                     "se.alpha","se.beta","se.sigma.alpha","se.sigma.beta")

