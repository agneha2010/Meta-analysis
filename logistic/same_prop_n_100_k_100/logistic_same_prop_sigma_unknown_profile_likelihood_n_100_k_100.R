library(locfit)
library(mvtnorm)
library(maxLik)
library(dplyr)

n=100
k=100
true.beta <- 0.5 ; true.alpha <- 0.5
true.sigma.beta <- 0.5 ; true.sigma.alpha <- 0.5
prop.trt <- 0.3
true.theta <- c(true.beta,true.sigma.beta,true.alpha,true.sigma.alpha)
x <- list() 
for (i in 1:k){
  x[[i]] <- rbinom(n,1,prop.trt)
}

comb.est <- function(n,k){
  #n=100; k=50
  beta.k <- rnorm(k,true.beta,sqrt(true.sigma.beta))
  alpha.k <- rnorm(k,true.alpha,sqrt(true.sigma.alpha))
  data <- list() 
  sigma.beta.hat.k.m <- c(0) ; beta.hat.k.m <- c(0)
  beta.hat.k <- list() ; sigma.beta.hat.k <-list()
  
  for (i in 1:k){
    #x[[i]] <- rbinom(n,1,prop.trt)
    p <- expit(alpha.k[i]*rep(1,n) + beta.k[i]*x[[i]])
    data[[i]] <- data.frame(x=x[[i]],y=rbinom(n,1,p))
    mylogit <- glm(y~x, data = data[[i]], family = "binomial")
    beta.hat.k[[i]] <- mylogit$coefficients
    sigma.beta.hat.k[[i]] <- vcov(mylogit)
    beta.hat.k.m[i] <- beta.hat.k[[i]][[2]]
    sigma.beta.hat.k.m[i] <- sigma.beta.hat.k[[i]][2,2]
  }
  
  comb.est.per.ipd <- function(per.ipd){
    #per.ipd <- 0.6
    k1 <- k*per.ipd
    k2 <- k*(1-per.ipd)
    beta.hat.k.ipd <- beta.hat.k[1:k1]
    sigma.beta.hat.k.ipd <- sigma.beta.hat.k[1:k1]
    beta.hat.k.meta <- beta.hat.k.m[(k1+1):k]
    sigma.beta.hat.k.meta <- sigma.beta.hat.k.m[(k1+1):k]
    logLikFun <- function(par){
      #par = c(0.1,2,0.5,2)
      alpha <- par[1]
      beta <- par[2]
      sigma.alpha <- par[3]
      sigma.beta <- par[4]
      
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
    constraints <- list(ineqA=matrix(c(0,0,1,0,0,0,0,1),2,4,byrow=T), ineqB=c(0,0))
    if(k1>0){
      mle <- maxLik(logLik=logLikFun,start=c(0.1,0.1,2,2),constraints=constraints,method="BFGS")
      }else{
      mle <- maxLik(logLik=logLikFun,start=c(alpha=true.alpha,beta=0.1,sigma.alpha=true.sigma.alpha,sigma.beta=2),
                    fixed=c("alpha","sigma.alpha"),constraints=constraints,method="BFGS")
      }
    gg <- summary(mle)
    est <- gg$estimate[,1]
    sd <- gg$estimate[,2]
    
    ####---- Using the closed form expressions
    theta.f = function(sigma.alpha,sigma.beta){
      #sigma.alpha= NA; sigma.beta = 0.488276
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
      if(k2==0){
        var = solve(Reduce("+",m))
        est = var%*%(Reduce("+",s1))
      } else if(k1==0){
        #tmp = 
        var = matrix(c(NA,NA,NA,solve(Reduce("+",l)[2,2])),nrow=2)
        est = c(NA,solve(Reduce("+",l)[2,2])%*%(Reduce("+",s2)[2]))
      }else{
        tmp1 = if (length(Reduce("+",m))==0) 0 else Reduce("+",m) + if (length(Reduce("+",l))==0) 0 else Reduce("+",l)
        tmp2 <- if (length(Reduce("+",s1))==0) 0 else Reduce("+",s1) + if (length(Reduce("+",s2))==0) 0 else Reduce("+",s2)
        var <- solve(tmp1)
        est <- var%*%tmp2
      }
      return(list(est,var))
    }
    
    logLikFun.1 <- function(par){
      #par = c(2,2)
      sigma.alpha <- par[1]
      sigma.beta <- par[2]
      theta = theta.f(sigma.alpha,sigma.beta)
      alpha = theta[[1]][1]
      beta = theta[[1]][2]
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
    constraints <- list(ineqA=matrix(c(1,0,0,1),2,2,byrow=T), ineqB=c(0,0))
    if(k1>0){
      mle.1 <- maxLik(logLik=logLikFun.1,start=c(2,2),constraints=constraints,method="BFGS")
    }else{
      mle.1 <- maxLik(logLik=logLikFun.1,start=c(sigma.alpha=true.sigma.alpha,sigma.beta=2),
                    fixed=c("sigma.alpha"),constraints=constraints,method="BFGS")
    }
    gg.1 <- summary(mle.1)
    est.1 <- gg.1$estimate[,1]
    sd.1 <- gg.1$estimate[,2]
    if(k1==0){
      est.1[1] = NA
      sd.1[1] = NA
    }
    est.2 = theta.f(est.1[1],est.1[2])
    est.final = c(est.2[[1]][1],est.2[[1]][2],est.1[1],est.1[2])
    sd.final = c(sqrt(est.2[[2]][1,1]),sqrt(est.2[[2]][2,2]),sd.1[1],sd.1[2])
    mm = cbind(c(est,sd),c(est.final,sd.final))
    return(mm)
  }
  output <- mapply(comb.est.per.ipd,prop.ipd)
  return(output)
}

nsim <- 500
prop.ipd <- c(0,0.2,0.4,0.6,0.8,1)
result <- list()
for(i in 1:nsim){
  result[[i]] <- comb.est(n,k)
  print(paste0("simulation ", i))
}
ipd0 <- matrix(0, nrow=nsim, ncol=8)
ipd20 <- matrix(0, nrow=nsim, ncol=8)
ipd40 <- matrix(0, nrow=nsim, ncol=8)
ipd60 <- matrix(0, nrow=nsim, ncol=8)
ipd80 <- matrix(0, nrow=nsim, ncol=8)
ipd100 <- matrix(0, nrow=nsim, ncol=8)
for(i in 1:nsim){
  ipd0[i,] <- result[[i]][1:8,1]
  ipd20[i,] <- result[[i]][1:8,2]
  ipd40[i,] <- result[[i]][1:8,3]
  ipd60[i,] <- result[[i]][1:8,4]
  ipd80[i,] <- result[[i]][1:8,5]
  ipd100[i,] <-result[[i]][1:8,6]
}
ipd0 <- as.data.frame(ipd0)
ipd20 <- as.data.frame(ipd20)
ipd40 <- as.data.frame(ipd40)
ipd60 <- as.data.frame(ipd60)
ipd80 <- as.data.frame(ipd80)
ipd100 <- as.data.frame(ipd100)
ipd0 <- ipd0 %>% 
  filter_all(all_vars(!is.infinite(.)))
ipd20 <- ipd20 %>% 
  filter_all(all_vars(!is.infinite(.)))
ipd40 <- ipd40 %>% 
  filter_all(all_vars(!is.infinite(.)))
ipd60 <- ipd60 %>% 
  filter_all(all_vars(!is.infinite(.)))
ipd80 <- ipd80 %>% 
  filter_all(all_vars(!is.infinite(.)))
ipd100 <- ipd100 %>% 
  filter_all(all_vars(!is.infinite(.)))
theta.hat.ipd0 <- c(colMeans(ipd0,na.rm=T)[1:4],colMeans(ipd0^2,na.rm=T)[5:8])
theta.hat.ipd20 <- c(colMeans(ipd20,na.rm=T)[1:4],colMeans(ipd20^2,na.rm=T)[5:8])
theta.hat.ipd40 <- c(colMeans(ipd40,na.rm=T)[1:4],colMeans(ipd40^2,na.rm=T)[5:8])
theta.hat.ipd60 <- c(colMeans(ipd60,na.rm=T)[1:4],colMeans(ipd60^2,na.rm=T)[5:8])
theta.hat.ipd80 <- c(colMeans(ipd80,na.rm=T)[1:4],colMeans(ipd80^2,na.rm=T)[5:8])
theta.hat.ipd100 <-c(colMeans(ipd100,na.rm=T)[1:4],colMeans(ipd100^2,na.rm=T)[5:8])
est <- t(data.frame(theta.hat.ipd0,theta.hat.ipd20, theta.hat.ipd40, theta.hat.ipd60,
                    theta.hat.ipd80, theta.hat.ipd100))
colnames(est)[1:4] <- c("alpha","beta","sigma.alpha","sigma.beta")
true.theta.mat <- matrix(kronecker(true.theta,rep(1,6)),byrow=F,nrow=6)
bias <- est[,1:4]-true.theta.mat
var <- est[,5:8]
sd <- rbind(apply(ipd0[,1:4],MARGIN=2,sd),
            apply(ipd20[,1:4],MARGIN=2,sd),
            apply(ipd40[,1:4],MARGIN=2,sd),
            apply(ipd60[,1:4],MARGIN=2,sd),
            apply(ipd80[,1:4],MARGIN=2,sd),
            apply(ipd100[,1:4],MARGIN=2,sd))

mse <- bias^2+sd^2
mse.com <- matrix(kronecker(mse[6,],rep(1,6)),byrow=F,nrow=6)
rel.eff <- mse.com/mse

####---- Using the closed forms
ipd0.1 <- matrix(0, nrow=nsim, ncol=8)
ipd20.1 <- matrix(0, nrow=nsim, ncol=8)
ipd40.1 <- matrix(0, nrow=nsim, ncol=8)
ipd60.1 <- matrix(0, nrow=nsim, ncol=8)
ipd80.1 <- matrix(0, nrow=nsim, ncol=8)
ipd100.1 <- matrix(0, nrow=nsim, ncol=8)
for(i in 1:nsim){
  ipd0.1[i,] <- result[[i]][9:16,1]
  ipd20.1[i,] <- result[[i]][9:16,2]
  ipd40.1[i,] <- result[[i]][9:16,3]
  ipd60.1[i,] <- result[[i]][9:16,4]
  ipd80.1[i,] <- result[[i]][9:16,5]
  ipd100.1[i,] <-result[[i]][9:16,6]
}
ipd0.1 <- as.data.frame(ipd0.1)
ipd20.1 <- as.data.frame(ipd20.1)
ipd40.1 <- as.data.frame(ipd40.1)
ipd60.1 <- as.data.frame(ipd60.1)
ipd80.1 <- as.data.frame(ipd80.1)
ipd100.1 <- as.data.frame(ipd100.1)
ipd0.1 <- ipd0.1 %>% 
  filter_all(all_vars(!is.infinite(.)))
ipd20.1 <- ipd20.1 %>% 
  filter_all(all_vars(!is.infinite(.)))
ipd40.1 <- ipd40.1 %>% 
  filter_all(all_vars(!is.infinite(.)))
ipd60.1 <- ipd60.1 %>% 
  filter_all(all_vars(!is.infinite(.)))
ipd80.1 <- ipd80.1 %>% 
  filter_all(all_vars(!is.infinite(.)))
ipd100.1 <- ipd100.1 %>% 
  filter_all(all_vars(!is.infinite(.)))
theta.hat.ipd0.1 <- c(colMeans(ipd0.1,na.rm=T)[1:4],colMeans(ipd0.1^2,na.rm=T)[5:8])
theta.hat.ipd20.1 <- c(colMeans(ipd20.1,na.rm=T)[1:4],colMeans(ipd20.1^2,na.rm=T)[5:8])
theta.hat.ipd40.1 <- c(colMeans(ipd40.1,na.rm=T)[1:4],colMeans(ipd40.1^2,na.rm=T)[5:8])
theta.hat.ipd60.1 <- c(colMeans(ipd60.1,na.rm=T)[1:4],colMeans(ipd60.1^2,na.rm=T)[5:8])
theta.hat.ipd80.1 <- c(colMeans(ipd80.1,na.rm=T)[1:4],colMeans(ipd80.1^2,na.rm=T)[5:8])
theta.hat.ipd100.1 <-c(colMeans(ipd100.1,na.rm=T)[1:4],colMeans(ipd100.1^2,na.rm=T)[5:8])
est.1 <- t(data.frame(theta.hat.ipd0.1,
                      theta.hat.ipd20.1, 
                      theta.hat.ipd40.1, 
                      theta.hat.ipd60.1,
                      theta.hat.ipd80.1, 
                      theta.hat.ipd100.1))
colnames(est.1)[1:4] <- c("alpha","beta","sigma.alpha","sigma.beta")
true.theta.mat <- matrix(kronecker(true.theta,rep(1,6)),byrow=F,nrow=6)
bias.1 <- est.1[,1:4]-true.theta.mat
var.1 <- est.1[,5:8]
sd.1 <- rbind(apply(ipd0.1[,1:4],MARGIN=2,sd),
              apply(ipd20.1[,1:4],MARGIN=2,sd),
              apply(ipd40.1[,1:4],MARGIN=2,sd),
              apply(ipd60.1[,1:4],MARGIN=2,sd),
              apply(ipd80.1[,1:4],MARGIN=2,sd),
              apply(ipd100.1[,1:4],MARGIN=2,sd))

mse.1 <- bias.1^2+sd.1^2
mse.com.1 <- matrix(kronecker(mse.1[6,],rep(1,6)),byrow=F,nrow=6)
rel.eff.1 <- mse.com.1/mse.1


#sum(ifelse(is.na(ipd20),1,0))-2
 # setwd("~/Desktop/Backup/meta_analysis/manuscript_material/logistic/unknown_sigma_profile_likelihood/simulation_n_100_k_100")
 # write.csv(ipd0, file = "ipd0.csv")
 # write.csv(ipd20, file = "ipd20.csv")
 # write.csv(ipd40, file = "ipd40.csv")
 # write.csv(ipd60, file = "ipd60.csv")
 # write.csv(ipd80, file = "ipd80.csv")
 # write.csv(ipd100, file = "ipd100.csv")
 # write.csv(est, file = "est.csv")
 # write.csv(bias, file = "bias.csv")
 # write.csv(var, file = "var_from_model.csv")
 # write.csv(mse, file = "mse.csv")
 # write.csv(rel.eff, file = "rel_eff.csv")
 # write.csv(sd, file = "sd.csv")
 # 
 # write.csv(ipd0.1, file = "ipd0.1.csv")
 # write.csv(ipd20.1, file = "ipd20.1.csv")
 # write.csv(ipd40.1, file = "ipd40.1.csv")
 # write.csv(ipd60.1, file = "ipd60.1.csv")
 # write.csv(ipd80.1, file = "ipd80.1.csv")
 # write.csv(ipd100.1, file = "ipd100.1.csv")
 # write.csv(est.1, file = "est.1.csv")
 # write.csv(bias.1, file = "bias.1.csv")
 # write.csv(var.1, file = "var_from_model.1.csv")
 # write.csv(mse.1, file = "mse.1.csv")
 # write.csv(rel.eff.1, file = "rel_eff.1.csv")
 # write.csv(sd.1, file = "sd.1.csv")
 # 
