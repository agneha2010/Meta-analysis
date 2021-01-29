library(openxlsx)
library(metafor)
library(lme4)
library(reshape2)
library(dplyr)
library(locfit)
library(mvtnorm)
library(maxLik)
library(dplyr)

mi.data = read.xlsx("~/Desktop/Backup/meta_analysis/manuscript_material/real_data/glm/data.xlsx")
k = nrow(mi.data)
mi.data$alive.control = mi.data$total.control-mi.data$death.control
mi.data$alive.treatment = mi.data$total.treatment-mi.data$death.treatment

#sum(dat.ma$yi/(dat.ma$vi+0.0103))/sum(1/(dat.ma$vi+0.0103)) # -0.2496771
#sqrt(1/sum(1/(dat.ma$vi+0.0103))) # 0.05812645

#rma.glmm(measure="OR", ai=mi.data$death.treatment,
#         bi=mi.data$alive.treatment,ci=mi.data$death.control,
#         di=mi.data$alive.control,
#         data=dat.ma, model="UM.RS", verbose=TRUE)

data = list()
sigma.beta.hat.k.m <- c(0) ; beta.hat.k.m <- c(0); alpha.hat.k.m = c(0)
beta.hat.k <- list() ; sigma.beta.hat.k <-list()

for (i in 1:k){
  #i=1
  x = c(rep(1,mi.data$death.control[i]),rep(2,mi.data$death.treatment[i]),
        rep(1,mi.data$alive.control[i]),rep(2,mi.data$alive.treatment[i]))
  y = c(rep(1,mi.data$death.control[i]+mi.data$death.treatment[i]),
        rep(0,mi.data$alive.control[i]+mi.data$alive.treatment[i]))
  data[[i]] = data.frame(x=x,y=y)
  mylogit <- glm(y~x, data = data[[i]], family = "binomial")
  beta.hat.k[[i]] <- mylogit$coefficients
  sigma.beta.hat.k[[i]] <- vcov(mylogit)
  alpha.hat.k.m[i] <- beta.hat.k[[i]][[1]]
  beta.hat.k.m[i] <- beta.hat.k[[i]][[2]]
  sigma.beta.hat.k.m[i] <- sigma.beta.hat.k[[i]][2,2]
}
ini.sigma.alpha = var(alpha.hat.k.m)
ini.sigma.beta = var(beta.hat.k.m)
ini.beta = mean(beta.hat.k.m)
ini.alpha = mean(alpha.hat.k.m)

comb.est.per.ipd <- function(per.ipd){
  #per.ipd <- 0.2
  k1 <- round(k*per.ipd,0)
  k2 <- round(k*(1-per.ipd),0)
  beta.hat.k.ipd <- beta.hat.k[1:k1]
  sigma.beta.hat.k.ipd <- sigma.beta.hat.k[1:k1]
  beta.hat.k.meta <- beta.hat.k.m[(k1+1):k]
  sigma.beta.hat.k.meta <- sigma.beta.hat.k.m[(k1+1):k]
  logLikFun <- function(par){
    #par = c(0.1,2,0.1,2)
    beta <- par[1]
    sigma.beta <- par[2]
    alpha <- par[3]
    sigma.alpha <- par[4]
    s1 <- c(0) ; s2 <- c(0)
    if(k1>0){
      for (i in 1:k1){
        m <- matrix(c(sigma.alpha,0,0,sigma.beta),nrow=2) + as.matrix(sigma.beta.hat.k.ipd[[i]])
        s1[i] <- dmvnorm(beta.hat.k.ipd[[i]], mean = c(alpha,beta), sigma = m, log = TRUE)
      }
    }
    if(k2>0){
      for (i in 1:k2){
        l <- sqrt(sigma.beta + sigma.beta.hat.k.meta[i])
        s2[i] <- dnorm(beta.hat.k.meta[i], mean = beta, sd = l, log = TRUE)
      }
    }
    return(sum(s1,na.rm = T) + sum(s2, na.rm = T))
  }
  constraints <- list(ineqA=matrix(c(0,1,0,0,0,0,0,1),2,4,byrow=T), ineqB=c(0,0))
  if(k1>0){
    mle <- maxLik(logLik=logLikFun,start=c(ini.beta,ini.sigma.beta,ini.alpha,ini.sigma.alpha),
                  constraints=constraints,method="BFGS")
    #mle <- maxLik(logLik=logLikFun,start=c(-0.24958595,0.01028180,-1.913599793,0.217346351),constraints=constraints,method="BFGS")
    }else{
    mle <- maxLik(logLik=logLikFun,start=c(beta=ini.beta,sigma.beta=ini.sigma.beta,alpha=ini.alpha,sigma.alpha=ini.sigma.alpha),
                  fixed=c("alpha","sigma.alpha"),constraints=constraints)
    }
  #test = optim(c(ini.beta,100,ini.alpha,100), logLikFun, NULL, method = "L-BFGS-B",
  #             lower = c(-Inf,0,-Inf,0), upper = rep(Inf, 4),hessian=TRUE)
  #cov.hessian = solve(test$hessian)
  gg <- summary(mle)
  est <- gg$estimate[,1]
  sd <- gg$estimate[,2]
  return(c(est,sd))
}

prop.ipd <- c(0,0.2,0.4,0.6,0.8,1)
round(k*prop.ipd,0)
output <- round(as.data.frame(mapply(comb.est.per.ipd,prop.ipd)),4)
colnames(output) = paste0("prop.ipd_",prop.ipd)
rownames(output)[5:8] = c("se.beta","se.sigma.beta","se.alpha","se.sigma.alpha")
#write.csv(output,"output.csv")
###------- Verifying the all.ipd and all.ma through random effects model

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
summary(mylogit)

###--- all MA estimate of \beta
mi.data$alive.control = mi.data$total.control-mi.data$death.control
mi.data$alive.treatment = mi.data$total.treatment-mi.data$death.treatment
dat.ma = data.frame(escalc(measure="OR",ai=mi.data$death.treatment,
                           bi=mi.data$alive.treatment,ci=mi.data$death.control,
                           di=mi.data$alive.control),center=1:k)   
#dat.ma1_4 = dat.ma[1:22,]
meta.RE = rma(yi,vi,data=dat.ma,method = "ML")
meta.RE
