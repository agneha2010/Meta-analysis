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

comb.est.per.ipd <- function(ipd.seq){
  ma.seq = which(!1:k %in% ipd.seq)
  k1 <- length(ipd.seq)
  k2 <- length(ma.seq)
  if(k1>0){
    beta.hat.k.ipd <- beta.hat.k[ipd.seq]
    sigma.beta.hat.k.ipd <- sigma.beta.hat.k[ipd.seq]
  }
  if(k2>0){
    beta.hat.k.meta <- beta.hat.k.m[ma.seq]
    sigma.beta.hat.k.meta <- sigma.beta.hat.k.m[ma.seq]
  }
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
    }else{
    mle <- maxLik(logLik=logLikFun,start=c(beta=ini.beta,sigma.beta=ini.sigma.beta,alpha=ini.alpha,sigma.alpha=ini.sigma.alpha),
                  fixed=c("alpha","sigma.alpha"),constraints=constraints)
    }
  gg <- summary(mle)
  est <- gg$estimate[,1]
  sd <- gg$estimate[,2]
  return(c(est,sd))
}

prop.ipd <- c(0,0.2,0.4,0.6,0.8,1)
ipd.seq.list = list(NULL,c(1,2,17,19),c(1:3,17:22),1:13,1:18,1:22)
output <- round(as.data.frame(mapply(comb.est.per.ipd,ipd.seq.list)),4)
colnames(output) = paste0("prop.ipd_",prop.ipd)
rownames(output)[5:8] = c("se.beta","se.sigma.beta","se.alpha","se.sigma.alpha")
#write.csv(output,"output_ipd_seq.csv")

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

# Generalized linear mixed model fit by maximum likelihood (Laplace
#                                                           Approximation) [glmerMod]
# Family: binomial  ( logit )
# Formula: cbind(death, size - death) ~ 1 + treatment + (1 + treatment |  
#                                                          center)
# Data: dat
# 
# AIC      BIC   logLik deviance df.resid 
# 327.7    336.6   -158.8    317.7       39 
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -1.45612 -0.40887  0.02525  0.44137  1.39542 
# 
# Random effects:
#   Groups Name        Variance Std.Dev. Corr 
# center (Intercept) 0.254613 0.5046        
# treatment2  0.009585 0.0979   -0.45
# Number of obs: 44, groups:  center, 22
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -2.20413    0.11749 -18.760  < 2e-16 ***
#   treatment2  -0.24681    0.05888  -4.192 2.77e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr)
# treatment2 -0.364

#######-------- Dinchen's exact result for glm data All IPD
mylogit1 = glmer(cbind(death,size-death) ~ 1 + treatment + (1 + treatment | center),
                data = dat, family = binomial(link="logit"), nAGQ = 0)
summary(mylogit1)

# Generalized linear mixed model fit by maximum likelihood (Adaptive
#                                                           Gauss-Hermite Quadrature, nAGQ = 0) [glmerMod]
# Family: binomial  ( logit )
# Formula: cbind(death, size - death) ~ 1 + treatment + (1 + treatment |  
#                                                          center)
# Data: dat
# 
# AIC      BIC   logLik deviance df.resid 
# 327.7    336.6   -158.8    317.7       39 
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -1.45908 -0.41182  0.02264  0.43727  1.39214 
# 
# Random effects:
#   Groups Name        Variance Std.Dev. Corr 
# center (Intercept) 0.254678 0.5047        
# treatment2  0.009565 0.0978   -0.45
# Number of obs: 44, groups:  center, 22
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -2.19293    0.11697  -18.75  < 2e-16 ***
#   treatment2  -0.24740    0.05687   -4.35 1.36e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr)
# treatment2 -0.357

###--- all MA estimate of \beta
mi.data$alive.control = mi.data$total.control-mi.data$death.control
mi.data$alive.treatment = mi.data$total.treatment-mi.data$death.treatment
dat.ma = data.frame(escalc(measure="OR",ai=mi.data$death.treatment,
                           bi=mi.data$alive.treatment,ci=mi.data$death.control,
                           di=mi.data$alive.control),center=1:k)   
#dat.ma1_4 = dat.ma[1:22,]
#saveRDS(dat.ma,"dat_ma.rds")
meta.RE = rma(yi,vi,data=dat.ma,method = "ML")
meta.RE
