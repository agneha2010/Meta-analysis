library(openxlsx)
library(metafor)
library(lme4)
library(reshape2)
library(dplyr)

mi.data = read.xlsx("~/Desktop/Backup/meta_analysis/manuscript_material/real_data/glm/data.xlsx")
k = nrow(mi.data)

###------- estimating the heterogeneity through random effects model

###--- all IPD estimate of \beta
death = data.frame(melt(mi.data,measure.vars=c("death.control","death.treatment"),
                        id="center",variable.name="trt.cntl",value.name="death"),
                   id=1:(2*k))
size = data.frame(melt(mi.data,measure.vars=c("total.control","total.treatment"),
                       id="center",variable.name="trt.cntl",value.name="size"),
                   id=1:(2*k))
dat = merge(death,size,by="id",all.x = T)
dat = dat %>% select(center.x,trt.cntl.x,death,size)
colnames(dat)[1:2] = c("center","treatment")
levels(dat$treatment) = c("1","2")
dat$center = as.factor(dat$center)
mylogit = glmer(cbind(death,size-death) ~ 1 + treatment + (1 + treatment | center),
                data = dat, family = binomial(link="logit"),nAGQ = 1)
#trtment = ranef(mylogit)
#ranef(mylogit, condVar = TRUE)
#colMeans(ranef(mylogit1)$treatment[])
#var(ranef(mylogit1)$treatment[])
#log(sum(mi.data$death.treatment)*sum(mi.data$alive.control)/(sum(mi.data$death.control)*sum(mi.data$alive.treatment)))

###--- all MA estimate of \beta
mi.data$alive.control = mi.data$total.control-mi.data$death.control
mi.data$alive.treatment = mi.data$total.treatment-mi.data$death.treatment
dat.ma = data.frame(escalc(measure="OR",ai=mi.data$death.treatment,
             bi=mi.data$alive.treatment,ci=mi.data$death.control,
             di=mi.data$alive.control),center=1:k)
meta.RE = rma(yi,vi,data=dat.ma,method = "ML")
meta.RE
#sum(dat.ma$yi/(dat.ma$vi+0.0103))/sum(1/(dat.ma$vi+0.0103)) # -0.2496771
#sqrt(1/sum(1/(dat.ma$vi+0.0103))) # 0.05812645

#rma.glmm(measure="OR", ai=mi.data$death.treatment,
#         bi=mi.data$alive.treatment,ci=mi.data$death.control,
#         di=mi.data$alive.control,
#         data=dat.ma, model="UM.RS", verbose=TRUE)

data = list()
sigma.meta <- c(0) ; beta.meta <- c(0)
beta.hat.k <- list() ; sigma.beta.hat.k <-list()

for (i in 1:k){
  i=1
  x = c(rep(1,mi.data$death.control[i]),rep(2,mi.data$death.treatment[i]),
        rep(1,mi.data$alive.control[i]),rep(2,mi.data$alive.treatment[i]))
  y = c(rep(1,mi.data$death.control[i]+mi.data$death.treatment[i]),
        rep(0,mi.data$alive.control[i]+mi.data$alive.treatment[i]))
  data[[i]] = data.frame(x=x,y=y)
  mylogit <- glm(y~x, data = data[[i]], family = "binomial")
  beta.hat.k[[i]] <- mylogit$coefficients
  sigma.beta.hat.k[[i]] <- vcov(mylogit)
  beta.meta[i] <- beta.hat.k[[i]][[2]]
  sigma.meta[i] <- sigma.beta.hat.k[[i]][2,2]
}

comb.est.per.ipd <- function(per.ipd){
  per.ipd <- 0
  k1 <- k*per.ipd
  k2 <- k*(1-per.ipd)
  beta.hat.k.ipd <- beta.hat.k[1:k1]
  sigma.beta.hat.k.ipd <- sigma.beta.hat.k[1:k1]
  beta.hat.k.meta <- beta.meta[(k1+1):k]
  sigma.beta.hat.k.meta <- sigma.meta[(k1+1):k]
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
  tmp1 <- Reduce("+",m)+Reduce("+",l)
  tmp2 <- Reduce("+",s1)+Reduce("+",s2)
  var <- solve(tmp1)
  est <- var%*%tmp2
  var1 <- c(var[1,1],var[2,2])
  return(c(est,var1))
}

prop.ipd <- c(0,0.2,0.4,0.6,0.8,1)
output <- mapply(comb.est.per.ipd,prop.ipd)


