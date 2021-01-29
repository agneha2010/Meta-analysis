library(CombMSC)
objfn<-function(y,x){
  s = sum(x)
  m = sum(y*x)/s
  return(sum(x*(y-m)^2))
}

k = 40

y = runif(k)
x = runif(k)
k1 = 6
pnt = proc.time()
M1 = subsets(k,k1)
v1 = matrix(0, dim(M1)[1],1)
for(i in 1:dim(M1)[1]){
  v1[i] =objfn(y[M1[i,]], x[M1[i,]])
}
oset = M1[which.max(v1),]


#k2 = 2
# M2 = subsets(k,k2)
# v2 = matrix(0, dim(M2)[1],1)
# for(i in 1:dim(M2)[1]){
#   v2[i] = objfn(y[M2[i,]], x[M2[i,]])  
# }
# iset = M2[which.max(v2),]


pnt = proc.time()

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



sort(iset)
sort(oset)


# 2 25 26 28 45

