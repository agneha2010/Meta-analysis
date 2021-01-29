library(locfit)
library(mvtnorm)
library(maxLik)
library(dplyr)
library(invgamma)
library(CombMSC)

n=10
k=30
true.beta <- 1.5 ; true.alpha <- 0.5
sigma.alpha <- 0.025
sigma = 2.5
#sigma = 1/rgamma(k,rate = 5,shape=2)
true.theta <- c(true.beta,true.alpha)
#prop.trt <- c(rep(0.1,4),rep(0.2,1),rep(0.8,1),rep(0.9,4))
b <- sigma.alpha/sigma
a <- 1+n*b

y = rep(n/(sigma*a),k)

rmixbeta<-function(){
  p = runif(1)
  if(p<=0.5){
    x = rbeta(1,2,9)
  }
  if(p>0.5){
    x = rbeta(1,9,2)
  }
  return(x)
}
prop.trt = round(as.numeric(purrr::rerun(k, rmixbeta())),1)

objfn<-function(x,y){
  s = sum(y)
  m = sum(x*y)/s
  return(sum(y*(x-m)^2))
}

ma.ipd.seq.expr <- function(ipd.seq,k1){
  #k1=2
  #ipd.seq = c(1,7)
  pi.tilda = sum(y[ipd.seq]*prop.trt[ipd.seq])/sum(y[ipd.seq])
  var = 1/(sum(y[ipd.seq]*(prop.trt[ipd.seq]-pi.tilda)^2)+sum(n*prop.trt*(1-prop.trt)/sigma))
  return(var)
  }

#comb.est <- function(n,k){

  optimal.seq <- function(k1){
    #k1=5
    M1 = subsets(k,k1)
    
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
      if ((k1-2) == 0) break
      ind = sset[which.max(y[sset]*(prop.trt[sset]- m)^2/(sw + y[sset]))]
      iset<- union(iset,ind)
      sset = setdiff(sset,ind)
      m = (sw*m + y[ind]*prop.trt[ind])/(sw + y[ind])
      sw = sw + y[ind]
    }
    iset = sort(iset)
    ratio = objfn(prop.trt[iset], y[iset])/objfn(prop.trt[oset], y[oset])
    mat = cbind(oset,iset,var[which.min(var)],ma.ipd.seq.expr(iset,k1),var[which.max(var)],ratio)
    return(mat)
  }
  optimal.seq(10)
#  return()
#}
ma.ipd.seq.expr(1:k,k)
ma.ipd.seq.expr(NULL,0)
#comb.est(n,k)
