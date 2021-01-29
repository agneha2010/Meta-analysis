library(locfit)
library(mvtnorm)
library(maxLik)
library(dplyr)
library(invgamma)
library(CombMSC)

n=10
k=10
true.beta <- 1.5 ; true.alpha <- 0.5
sigma.alpha <- 0.025; sigma = 2.5
true.theta <- c(true.beta,true.alpha)
prop.trt <- c(rep(0.1,4),rep(0.2,1),rep(0.8,1),rep(0.9,4))
b <- sigma.alpha/sigma
a <- 1+n*b
#rgamma(k,scale = 5,shape=2)
y = rep(n/(sigma*a),k)

x <- list() 
for (i in 1:k){
  x[[i]] <- rbinom(n,1,prop.trt)
}

objfn<-function(x,y){
  s = sum(y)
  m = sum(x*y)/s
  return(sum(y*(x-m)^2))
}

ma.ipd.seq.expr <- function(ipd.seq,k1){   
  #ipd.seq = 1:k
  tot.seq <- 1:k
  k2 <- k-k1
  ma.seq <- subset(tot.seq,!(tot.seq)%in% ipd.seq)
  n.prop.meta <- cbind(prop.trt[ma.seq],1-prop.trt[ma.seq])
  n.prop.ipd <- cbind(prop.trt[ipd.seq],1-prop.trt[ipd.seq])
  var.meta <- n*sum(n.prop.meta[,1]*n.prop.meta[,2])
  var.ipd <- ((n^2)*b/a)*sum(n.prop.ipd[,1]*n.prop.ipd[,2])+
    (n/(a*k1))*(sum(n.prop.ipd[,1])*(k1-sum(n.prop.ipd[,1])))
  #####----- check var expression
  #var = sigma/(n*sum((n.prop.ipd[,1]-mean(n.prop.ipd[,1]))^2)/a+n*sum(prop.trt*(1-prop.trt)))
  #sigma/(var.meta+var.ipd)
  #####----------------------------
  if(k1==0){
    return(sigma/var.meta)}
  else if (k2==0){
    return(sigma/var.ipd)}
  else{
    return(sigma/(var.meta+var.ipd))}
}

#comb.est <- function(n,k){

  optimal.seq <- function(k1){
    #k1=2
    M1 = subsets(k,k1)
    # v1 = matrix(0, dim(M1)[1],1)
    # for(i in 1:dim(M1)[1]){
    #   v1[i] =objfn(prop.trt[M1[i,]], y[M1[i,]])
    # }
    # oset = sort(M1[which.max(v1),])
    
    ipd.seq.list = as.list(data.frame(t(M1)))
    var = mapply(ma.ipd.seq.expr,ipd.seq.list,k1=k1)
    oset = sort(M1[which.min(var),])
    
    v = rep(1,k)
    z = 1/y
    d = (prop.trt%*%t(v) - v%*%t(prop.trt))^2/(v%*%t(z) + z%*%t(v))
    iset =c(which(d ==max(d), arr.ind = T)[1,])
    sset = setdiff(seq(1,k), iset)
    sw = sum(y[iset])
    m = sum(prop.trt[iset]*y[iset])/sw
    #v = objfn(prop.trt[iset], y[iset])
    
    for(l in 1:(k1-2)){
      if (k1-2 == 0) break
      ind = sset[which.max(y[sset]*(prop.trt[sset]- m)^2/(sw + y[sset]))]
      iset<- union(iset,ind)
      sset = setdiff(sset,ind)
      m = (sw*m + y[ind]*prop.trt[ind])/(sw + y[ind])
      sw = sw + y[ind]
    }
    iset = sort(iset)
    ratio = objfn(prop.trt[iset], y[iset])/objfn(prop.trt[oset], y[oset])
    mat = cbind(oset,iset)
    return(list(var,mat,ratio))
  }
  optimal.seq(2)
#  return()
#}
ma.ipd.seq.expr(1:k,k)
ma.ipd.seq.expr(NULL,0)
#comb.est(n,k)
