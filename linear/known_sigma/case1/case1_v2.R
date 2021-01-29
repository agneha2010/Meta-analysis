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
sigma = 2.5
true.theta <- c(true.beta,true.alpha)
b <- sigma.alpha/sigma
a <- 1+n*b

y = n/(sigma*a)

p <- c(0.1, 0.2, 0.3, 0.3, 0.3, 0.5, 0.6, 0.6, 0.8, 0.8)

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

comb.est <- function(prop.trt.sim){
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
prop.trt.sim = matrix(NA,nrow=nsim,ncol = k)
for(i in 1:nsim){
  print(paste0("simulation ", i))
  set.seed(i)
  prop.trt.sim[i,] = round(as.numeric(purrr::rerun(k, rmixbeta())),1)
  result[[i]] <- round(comb.est(prop.trt.sim[i,]),4)
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
    t[j]=ifelse(sum(oset1 %in% iset1)==(j-1),1,0)
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

### number of places in which they match == 0,1,...,10
colnames(ipd2) = c(paste0("match_count_",0:2),"min_var_oset","var_iset",
                   "max_var_oset","ratio_obj_iset_oset")
colnames(ipd3) = c(paste0("match_count_",0:3),"min_var_oset","var_iset",
                   "max_var_oset","ratio_obj_iset_oset")
colnames(ipd4) = c(paste0("match_count_",0:4),"min_var_oset","var_iset",
                   "max_var_oset","ratio_obj_iset_oset")
colnames(ipd5) = c(paste0("match_count_",0:5),"min_var_oset","var_iset",
                   "max_var_oset","ratio_obj_iset_oset")
colnames(ipd6) = c(paste0("match_count_",0:6),"min_var_oset","var_iset",
                   "max_var_oset","ratio_obj_iset_oset")
colnames(ipd7) = c(paste0("match_count_",0:7),"min_var_oset","var_iset",
                   "max_var_oset","ratio_obj_iset_oset")
colnames(ipd8) = c(paste0("match_count_",0:8),"min_var_oset","var_iset",
                   "max_var_oset","ratio_obj_iset_oset")
colnames(ipd9) = c(paste0("match_count_",0:9),"min_var_oset","var_iset",
                   "max_var_oset","ratio_obj_iset_oset")
colnames(ipd10) = c(paste0("match_count_",0:10),"min_var_oset","var_iset",
                   "max_var_oset","ratio_obj_iset_oset")
count_matches = cbind(colSums(ipd10[,1:11]),colSums(ipd9[,1:10]),colSums(ipd8[,1:9]),
                      colSums(ipd7[,1:8]),colSums(ipd6[,1:7]),colSums(ipd5[,1:6]),
                      colSums(ipd4[,1:5]),colSums(ipd3[,1:4]),colSums(ipd2[,1:3]))
colnames(count_matches) = paste0("k1_",10:2)

mean_ratio = c(mean(ipd10[,15]),mean(ipd9[,14]),mean(ipd8[,13]),
          mean(ipd7[,12]),mean(ipd6[,11]),mean(ipd5[,10]),
          mean(ipd4[,9]),mean(ipd3[,8]),mean(ipd2[,7]))
count_matches = rbind(count_matches,mean_ratio)

#write.csv(count_matches,"count_matches.csv")
#write.csv(prop.trt.sim,"prop.trt.sim.csv")
# write.csv(ipd2,"ipd2.csv")
# write.csv(ipd3,"ipd3.csv")
# write.csv(ipd4,"ipd4.csv")
# write.csv(ipd5,"ipd5.csv")
# write.csv(ipd6,"ipd6.csv")
# write.csv(ipd7,"ipd7.csv")
# write.csv(ipd8,"ipd8.csv")
# write.csv(ipd9,"ipd9.csv")
# write.csv(ipd10,"ipd10.csv")


