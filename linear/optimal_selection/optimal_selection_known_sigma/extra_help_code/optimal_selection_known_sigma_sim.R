library(CombMSC)
objfn<-function(y,x){
  s = sum(x)
  m = sum(y*x)/s
  return(sum(x*(y-m)^2))
}

pnt = proc.time()

nsim = 500
k = 50
k1 = 5
M1 = subsets(k,k1)
x.data = matrix(0,nrow=nsim,ncol=k)
y.data = matrix(0,nrow=nsim,ncol=k)
oset.store = matrix(0,nrow=nsim,ncol=k1)
iset.store = matrix(0,nrow=nsim,ncol=k1)
ratio = c()

for (j in 1:nsim){
  #j=1
  print(j)
  y = runif(k)
  x = runif(k)
  x.data[j,] = x
  y.data[j,] = y
  
  v1 = matrix(0, dim(M1)[1],1)
  for(i in 1:dim(M1)[1]){
    v1[i] =objfn(y[M1[i,]], x[M1[i,]])
  }
  oset = M1[which.max(v1),]
  
  v = rep(1,k)
  z = 1/x
  a = (y%*%t(v) - v%*%t(y))^2/(v%*%t(z) + z%*%t(v))
  iset =c(which(a ==max(a), arr.ind = T)[1,])
  sset = setdiff(seq(1,k), iset)
  sw = sum(x[iset])
  m = sum(y[iset]*x[iset])/sw
  v = objfn(y[iset], x[iset])
  
  for(l in 1:(k1-2)){
    ind = sset[which.max(x[sset]*(y[sset]- m)^2/(sw + x[sset]))]
    iset<- union(iset,ind)
    sset = setdiff(sset,ind)
    m = (sw*m + x[ind]*y[ind])/(sw + x[ind])
    sw = sw + x[ind]
  }
  oset.store[j,] = sort(oset)
  iset.store[j,] = sort(iset)
  ratio[j] = objfn(y[iset], x[iset])/objfn(y[oset], x[oset])
}

proc.time() - pnt
#user   system  elapsed 
#1612.088    3.346 1617.246
t = matrix(NA,nrow=nsim,ncol=(k1+1))
for(i in 1:nsim){
   for(j in 1:(k1+1)){
    t[i,j]=ifelse(sum(oset.store[i,]==iset.store[i,])==(j-1),1,0)
   }
}
per = colSums(t)/nsim*100 # 0.93
# 0.002 0.012 0.012 0.018 0.026 0.930
store = cbind(oset.store,round(ratio,3),iset.store)
colnames(store) = c(paste("optim",1:5),"ratio",paste("sub optim",1:5))
write.csv(store,"seq_ratio_store.csv")
write.csv(x.data,"x_data.csv")
write.csv(y.data,"y_data.csv")


