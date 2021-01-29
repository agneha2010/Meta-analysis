library(dplyr)
library(lme4)
#library(arm)
library(CombMSC)
library(reshape2)

data = readRDS("~/Desktop/Backup/meta_analysis/manuscript_material/real_data/lm/F00007762-WV6_Data_R_v20180912.rds")
# dim 89565 440
cntry.freq = as.data.frame(ftable(data$V2))
length(unique(data$V2))
cntry.regn.freq = as.data.frame(ftable(data$V2A))
year.freq = as.data.frame(ftable(data$V1))
life.satisfaction.freq = as.data.frame(ftable(data$V23))
financial.satisfaction.freq = as.data.frame(ftable(data$V55))
#V23 Satisfaction with your life
#V59 Satisfaction with financial situation of household
#V55 How much freedom of choice and control over own life
#V240 Sex
#V242 Age
mydata = data %>% filter(!V59 %in% c(-5,-2,-1)) %>%
                  filter(!V23 %in% c(-5,-2,-1)) %>% 
                  filter(!V240 %in% c(-5,-2,-1)) %>% 
                  filter(!V242 %in% c(-5,-3,-2,-1)) %>% 
                  filter(!V254 %in% c(-5,-4,-2,-1)) %>% 
                  filter(!V55 %in% c(-5,-2,-1)) %>% 
                  mutate(x = ifelse(V59<=5,0,1),age = case_when(V242<=20~"<=20",
                                                                V242>20 & V242<=40 ~"20-40",
                                                                V242>40 & V242<=60 ~"40-60",
                                                                V242>60 & V242<=80 ~"60-80",
                                                                V242>80~">80"),
                         x1 = ifelse(V55<=5,0,1)) %>% 
                  select(V2A,V23,x,x1,V240,age,V254) %>%
                  rename(y=V23,region=V2A,sex=V240,ethnicity=V254) 

mydata$region.sex.age.ethnicity = paste(paste0(mydata$region),"_",paste0(mydata$sex),
                              "_",paste0(mydata$age),"_",paste0(mydata$ethnicity))

### subsetting data for USA
mydata = mydata[mydata$region==840,]
code = unique(mydata$region.sex.age.ethnicity)
k = length(code)
data.list = list()
n = c()
prop.trt.mat = matrix(NA,ncol=2,nrow=k)
for(i in 1:k){
  #i=1
  data.list[[i]] = mydata[mydata$region.sex.age.ethnicity==code[i],]
  n[i] = dim(data.list[[i]])[1]
  prop.trt.mat[i,] = prop.table(table(data.list[[i]][,3]))
}
#write.csv(prop.trt.mat,"prop_trt_mat.csv")
prop.trt = prop.trt.mat[,1]
prop.trt = prop.trt[-which(prop.trt==1)]
mydata$region.sex.age.ethnicity = as.factor(mydata$region.sex.age.ethnicity)
mylm = lm(y~x,data = mydata)
summary(mylm)
mylmer = lmer(y ~ x + (1| region.sex.age.ethnicity),data = mydata)
summary(mylmer)
#ranef(mylmer)  fixef(mylmer) confint(mylmer) print(mylmer)

sigma.alpha = 0.04027
sigma = 2.72339

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
test = optimal.seq(5) # 0.005626499 
ma.ipd.seq.expr(1:k,k) # 0.00552328
ma.ipd.seq.expr(NULL,0) # 0.005705427
