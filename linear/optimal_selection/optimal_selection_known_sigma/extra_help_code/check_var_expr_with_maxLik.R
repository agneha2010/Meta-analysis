library(locfit)
library(mvtnorm)
library(maxLik)
library(dplyr)

n=20
k=100
beta <- 1.5 ;alpha <- 0.5; sigma <- 2.5
sigma.alpha <- 0.025
true.theta <- c(true.beta,true.alpha)
#set.seed(100)
prop.trt <- c(rep(0.1,40),rep(0.2,10),rep(0.8,10),rep(0.9,40))

x <- list() 
for (i in 1:k){
  x[[i]] <- rbinom(n,1,prop.trt)
}

#comb.est <- function(n,k){
  #n=1000; k=50
  alpha.k <- rnorm(k,true.alpha,sqrt(true.sigma.alpha))
  data <- list() 
  sigma.beta.hat.k.m <- c() ; beta.hat.k.m <- c()
  beta.hat.k <- list() ; sigma.beta.hat.k <-list()
  n.prop <- matrix(0,nrow=k,ncol=2)
  n.prop[,1] <- prop.trt
  n.prop[,2] <- 1-prop.trt
  sigma.alpha <- true.sigma.alpha
  sigma <- true.sigma
  
  for (i in 1:k){
    #i=4
    #set.seed(1000)
    y <- alpha.k[i]+true.beta*x[[i]]+rnorm(n,mean=0,sd=sqrt(true.sigma))
    data[[i]] <- data.frame(x=x[[i]],y=y)
    mylm <- lm(y~x, data = data[[i]])
    beta.hat.k.m[i] <- mylm$coefficients[2]
    sigma.beta.hat.k.m[i] <- vcov(mylm)[2,2]
  }
  ma.ipd.seq.expr <- function(ma.seq){
    #ma.seq = ma.seq.list20$X1
    ipd.seq <- subset(tot.seq,!(tot.seq)%in% ma.seq)
    k1 <- length(ipd.seq)
    k2 <- length(ma.seq)
    n.prop.meta <- n.prop[ma.seq,]
    n.prop.ipd <- n.prop[ipd.seq,]
    var.meta <- n*sum(n.prop.meta[,1]*n.prop.meta[,2])
    b <- sigma.alpha/sigma
    a <- 1+n*b
    y = rep(n/(sigma*a),k)
    pi.tilda = sum(y[ipd.seq]*prop.trt[ipd.seq])/sum(y[ipd.seq])
    var = 1/(sum(y[ipd.seq]*(prop.trt[ipd.seq]-pi.tilda)^2)+sum(n*prop.trt*(1-prop.trt)/sigma))
    return(var)
    var.ipd <- ((n^2)*b/a)*sum(n.prop.ipd[,1]*n.prop.ipd[,2])+
      (n/(a*k1))*(sum(n.prop.ipd[,1])*sum(n.prop.ipd[,2]))
    if(k1==0){
      return(sigma/var.meta)}
    else if (k2==0){
      return(sigma/var.ipd)}
    else{
      return(sigma/(var.meta+var.ipd))}
  }
  
  ma.ipd.seq.maxlik <- function(ma.seq){
    #ma.seq = ma.seq.list80$X1
    ipd.seq <- subset(tot.seq,!(tot.seq)%in% ma.seq)
    k1 <- length(ipd.seq)
    k2 <- length(ma.seq)
    data.ipd <- data[ipd.seq]
    if(k2>0){
      beta.hat.k.meta = beta.hat.k.m[ma.seq]
      sigma.beta.hat.k.meta <- sigma.beta.hat.k.m[ma.seq]
    }

    logLikFun <- function(par){
      #par = c(0.1,0.1)
      beta <- par[1]
      alpha <- par[2]
      s1 <- c(0) ; s2 <- c(0)
      H <- sigma.alpha*matrix(1,nrow=n,ncol=n)+sigma*diag(1,n)
      if(k1>0){
        for (i in 1:k1){
          m <- cbind(rep(1,n),data.ipd[[i]][[1]])%*%c(alpha,beta)
          s1[i] <- dmvnorm(data.ipd[[i]][[2]], mean = m, sigma = H, log = TRUE)
        }
      }
      if(k2>0){
        for (i in 1:k2){
          l <- sqrt(sigma.beta.hat.k.meta[i])
          s2[i] <- dnorm(beta.hat.k.meta[i], mean = beta, sd = l, log = TRUE)
        }
      }
      return(sum(s1,na.rm = T) + sum(s2, na.rm = T))
    }
    #constraints = list(ineqA=matrix(c(0,0,1,0,0,0,0,1),2,4,byrow=T), ineqB=c(0,0))
    if(k1>0){
      mle <- maxLik(logLik=logLikFun,start=c(0.1,0.1))
    }else{
      mle <- maxLik(logLik=logLikFun,start=c(beta=0.1,alpha=true.alpha),
                    fixed=c("alpha"))
    }
    gg <- summary(mle)
    est <- gg$estimate[,1]
    var <- (gg$estimate[,2])^2
    return(c(est,var))
  }
#  output <- mapply(comb.est.per.ipd,prop.ipd)
#  return(output)
  tot.seq <- 1:k
  ma.seq.list100 = as.list(data.frame(combn(k,k)))
  ma.seq.list80 = as.list(data.frame(combn(k,8)))
  ma.seq.list60 = as.list(data.frame(combn(k,6)))
  ma.seq.list40 = as.list(data.frame(combn(k,4)))
  ma.seq.list20 = as.list(data.frame(combn(k,2)))
  ma.seq.list0 = as.list(data.frame(combn(k,0)))
  
  ipd0.maxlik = mapply(ma.ipd.seq.maxlik,ma.seq.list100)
  ipd20.maxlik = mapply(ma.ipd.seq.maxlik,ma.seq.list80)
  ipd40.maxlik = mapply(ma.ipd.seq.maxlik,ma.seq.list60)
  ipd60.maxlik = mapply(ma.ipd.seq.maxlik,ma.seq.list40)
  ipd80.maxlik = mapply(ma.ipd.seq.maxlik,ma.seq.list20)
  ipd100.maxlik = mapply(ma.ipd.seq.maxlik,ma.seq.list0)
  
  ipd0.expr = mapply(ma.ipd.seq.expr,ma.seq.list100)
  ipd20.expr = mapply(ma.ipd.seq.expr,ma.seq.list80)
  ipd40.expr = mapply(ma.ipd.seq.expr,ma.seq.list60)
  ipd60.expr = mapply(ma.ipd.seq.expr,ma.seq.list40)
  ipd80.expr = mapply(ma.ipd.seq.expr,ma.seq.list20)
  ipd100.expr = mapply(ma.ipd.seq.expr,ma.seq.list0)
#  return(list(ipd0,ipd20,ipd40,ipd60,ipd80,ipd100))
#}

