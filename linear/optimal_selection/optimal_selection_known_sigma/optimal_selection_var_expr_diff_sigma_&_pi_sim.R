library(locfit)
library(mvtnorm)
library(maxLik)
library(dplyr)
library(invgamma)
library(CombMSC)

n = 10
k = 30
true.beta = 1.5 ; true.alpha = 0.5
sigma.alpha = 0.025
#sigma = 2.5
set.seed(1000)
sigma = 1/rgamma(k,rate = 5,shape=2)
true.theta <- c(true.beta,true.alpha)
b <- sigma.alpha/sigma
a <- 1+n*b

y = n/(sigma*a)

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

objfn<-function(x,y){
  s = sum(y)
  m = sum(x*y)/s
  return(sum(y*(x-m)^2))
}
ma.ipd.seq.expr <- function(ipd.seq,prop.trt,k1){
  #k1=k
  #ipd.seq = 1:k
  #prop.trt = prop.trt.sim
  pi.tilda = sum(y[ipd.seq]*prop.trt[ipd.seq])/sum(y[ipd.seq])
  var = 1/(sum(y[ipd.seq]*(prop.trt[ipd.seq]-pi.tilda)^2)+sum(n*prop.trt*(1-prop.trt)/sigma))
  return(var)
}

optimal.seq <- function(k1,subset.matrix,prop.trt){
  #k1=2
  #M1 = subsets(k,k1)
  M1 = subset.matrix
  ipd.seq.list = as.list(data.frame(t(M1)))
  var = mapply(ma.ipd.seq.expr,ipd.seq.list,MoreArgs = list(prop.trt=prop.trt,k1=k1))
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
  ratio = objfn(prop.trt[iset], y[iset])/objfn(prop.trt[oset], y[oset])
  mat = cbind(oset,sort(iset),var[which.min(var)],ma.ipd.seq.expr(iset,prop.trt,k1),var[which.max(var)],ratio)
  return(mat)
}

m2 = subsets(k,2)
m3 = subsets(k,3)
m4 = subsets(k,4)
m5 = subsets(k,5)
m6 = subsets(k,6)
m7 = subsets(k,7)
m8 = subsets(k,8)
m9 = subsets(k,9)
m10 = subsets(k,10)

comb.est <- function(k){
  prop.trt.sim = round(as.numeric(purrr::rerun(k, rmixbeta())),1)
  store = rbind(optimal.seq(2,m2,prop.trt = prop.trt.sim),optimal.seq(3,m3,prop.trt = prop.trt.sim),
                optimal.seq(4,m4,prop.trt = prop.trt.sim),optimal.seq(5,m5,prop.trt = prop.trt.sim),
                optimal.seq(6,m6,prop.trt = prop.trt.sim),optimal.seq(7,m7,prop.trt = prop.trt.sim),
                optimal.seq(8,m8,prop.trt = prop.trt.sim),optimal.seq(9,m9,prop.trt = prop.trt.sim),
                optimal.seq(10,m10,prop.trt = prop.trt.sim),ma.ipd.seq.expr(NULL,prop.trt.sim,0),
                ma.ipd.seq.expr(1:k,prop.trt.sim,k))
  return(store)
}
nsim = 100
result <- list()
for(i in 1:nsim){
  print(paste0("simulation ", i))
  result[[i]] <- round(comb.est(k),4)
}

ipd2 <- matrix(0, nrow=nsim, ncol=7)
ipd3 <- matrix(0, nrow=nsim, ncol=8)
ipd4 <- matrix(0, nrow=nsim, ncol=9)
ipd5 <- matrix(0, nrow=nsim, ncol=10)
ipd6 <- matrix(0, nrow=nsim, ncol=11)
ipd7 <- matrix(0, nrow=nsim, ncol=12)
ipd8 <- matrix(0, nrow=nsim, ncol=13)
ipd9 <- matrix(0, nrow=nsim, ncol=14)
ipd10 <- matrix(0, nrow=nsim, ncol=15)
ipd0 = c()
ipdk = c()

check.matches = function(oset1,iset1){
  #oset1 = result[[i]][1:2,1]
  #iset1 = result[[i]][1:2,2]
  k1 = length(oset1)
  t = c()
  for(j in 1:(k1+1)){
    t[j]=ifelse(sum(oset1==iset1)==(j-1),1,0)
  }
  return(t)
}

for(i in 1:nsim){
  ipd2[i,] <- c(check.matches(result[[i]][1:2,1],result[[i]][1:2,2]),result[[i]][1,3:6])
  ipd3[i,] <- c(check.matches(result[[i]][3:5,1],result[[i]][3:5,2]),result[[i]][3,3:6])
  ipd4[i,] <- c(check.matches(result[[i]][6:9,1],result[[i]][6:9,2]),result[[i]][6,3:6])
  ipd5[i,] <- c(check.matches(result[[i]][10:14,1],result[[i]][10:14,2]),result[[i]][10,3:6])
  ipd6[i,] <- c(check.matches(result[[i]][15:20,1],result[[i]][15:20,2]),result[[i]][15,3:6])
  ipd7[i,] <- c(check.matches(result[[i]][21:27,1],result[[i]][21:27,2]),result[[i]][21,3:6])
  ipd8[i,] <- c(check.matches(result[[i]][28:35,1],result[[i]][28:35,2]),result[[i]][28,3:6])
  ipd9[i,] <- c(check.matches(result[[i]][36:44,1],result[[i]][36:44,2]),result[[i]][38,3:6])
  ipd10[i,] <- c(check.matches(result[[i]][45:54,1],result[[i]][45:54,2]),result[[i]][45,3:6])
  ipd0[i] <- result[[i]][55,1]
  ipdk[i] <- result[[i]][56,1]
}

