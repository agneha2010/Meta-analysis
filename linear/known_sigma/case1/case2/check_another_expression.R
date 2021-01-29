
library(mvtnorm)
library(dplyr)
library(Matrix)

k <- 10
n <- 10
beta <- 1.5 ;alpha <- 0.5; sigma <- 2.5
sigma.alpha <- 0.025
true.theta <- c(alpha,beta)
set.seed(100)
alpha.k <- rnorm(k,alpha,sqrt(sigma.alpha))
b <- sigma.alpha/sigma

p <- c(0.1, 0.1, 0.1, 0.1, 0.2, 0.8, 0.9, 0.9, 0.9, 0.9)
p <- c(0.1, 0.2, 0.3, 0.3, 0.3, 0.5, 0.6, 0.6, 0.8, 0.8)

comb.est <- function(n,k){
  #n=50; k=10
  n.prop <- matrix(0,nrow=k,ncol=2)
  n.prop[,1] <- p
  n.prop[,2] <- 1-p
  
  ma.ipd.seq <- function(ma.seq){
    #ma.seq <- 1:k
    ipd.seq <- subset(tot.seq,!(tot.seq)%in% ma.seq)
    k1 <- length(ipd.seq)
    k2 <- length(ma.seq)
    n.prop.meta <- n.prop[ma.seq,]
    n.prop.ipd <- n.prop[ipd.seq,]
    var.ipd <- sum((n.prop.ipd[,1])^2) - (sum(n.prop.ipd[,1]))^2/k1
      return(var.ipd)
  }
  output <- mapply(ma.ipd.seq,ma.seq.list)
  return(output)
}

tot.seq <- 1:k
k2 <- 4
ma.seq.list <- as.list(data.frame(combn(k,k2)))
no.comb <- ncol(combn(k,k2))
#ma.seq.list <- c(ma.seq.list,list(all.meta=tot.seq),list(all.ipd=NULL))
mse <- comb.est(n,k)
rel.eff <- 0.1056068/mse
m <- data.frame(as.character(ma.seq.list),round(mse,4),round(rel.eff,4))
colnames(m) <- c("ma_seq","mse","rel.eff")

 order <- order(m[,2])
 m_ordered <- cbind(order,m[order,])
 m_ordered <- noquote(m_ordered)
 m_ordered
