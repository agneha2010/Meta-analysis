library(dplyr)
library(lme4)
#library(arm)
library(CombMSC)
#library(reshape2)

data = readRDS("~/Desktop/Backup/meta_analysis/manuscript_material/real_data/lm/F00007762-WV6_Data_R_v20180912.rds")
# dim 89565 440
data1 = data[data$V2A==840,]
#V23 Satisfaction with your life
#V59 Satisfaction with financial situation of household
#V55 How much freedom of choice and control over own life
#V240 Sex
#V242 Age
#V57 marital_status
#V254 Ethnicity
#V2A region
#V248_CS Education
mydata = data1 %>% filter(!V23 %in% c(-5,-4,-3,-2,-1)) %>% 
  filter(!V55 %in% c(-5,-4,-3,-2,-1)) %>% 
  #filter(V242 >= 40) %>%
  mutate(x = ifelse(V55<=5,0,1)) %>%
  rename(y=V23) %>%
  select(y,x) 
  
n_total =dim(mydata)[1]

cuts = floor(n_total/30)
temp = round(30/n_total,4)
prob = c(rep(temp,cuts-1),1-temp*(cuts-1))
#set.seed(1000)
idx <- sample(seq(1, cuts), size = nrow(mydata), replace = TRUE, prob = prob)
z = list()
for (i in 1:cuts)
  z[[i]] <- mydata[idx == i,]
mydata$idx = idx
k = cuts
n = c()
prop.trt.mat = matrix(NA,ncol=2,nrow=k)
for(i in 1:k){
  #i=1
  n[i] = dim(z[[i]])[1]
  prop.trt.mat[i,] = prop.table(table(z[[i]][,2]))
}
prop.trt = prop.trt.mat[,1]
inx = which(prop.trt==1)

n = n[-inx]
prop.trt = prop.trt[-inx]
k = length(prop.trt)
mydata = subset(mydata,!idx %in% inx)
store = data.frame(n=n,prop.trt=prop.trt)

mylm = lm(y~x,data = mydata)
mylmer = lmer(y ~ x + (1| idx),data = mydata)
summary(mylmer)
#ranef(mylmer)  fixef(mylmer) confint(mylmer) print(mylmer)

sigma.alpha = 0.006893    
sigma = 2.826807  

b <- sigma.alpha/sigma
a <- 1+n*b
y = n/(sigma*a)

objfn<-function(x,y){
  s = sum(y)
  m = sum(x*y)/s
  return(sum(y*(x-m)^2))
}

ma.ipd.seq.expr <- function(ipd.seq,k1){
  tot.seq <- 1:k
  k2 <- k-k1
  pi.tilda = sum(y[ipd.seq]*prop.trt[ipd.seq])/sum(y[ipd.seq])
  var = 1/(sum(y[ipd.seq]*(prop.trt[ipd.seq]-pi.tilda)^2)+sum(n*prop.trt*(1-prop.trt)/sigma))
  return(var)
}

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
  mat = cbind(oset,iset)
  return(list(var,mat,var[which.min(var)],ratio))
}
ma.ipd.seq.expr(1:k,k) # 0.0154988
ma.ipd.seq.expr(NULL,0) # 0.0208022
(ma.ipd.seq.expr(NULL,0)-ma.ipd.seq.expr(1:k,k))/ma.ipd.seq.expr(1:k,k)

test = optimal.seq(5) # 0.005626499 
optimal.seq(10)

