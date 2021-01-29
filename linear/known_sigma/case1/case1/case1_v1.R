
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

# p <- c(0.13, 0.15, 0.34, 0.42, 0.46, 0.55, 0.55, 0.56, 0.62, 0.76)
# p <- sort(round(runif(k,min=0,max=1),1))
p <- c(0.1, 0.2, 0.3, 0.3, 0.3, 0.5, 0.6, 0.6, 0.8, 0.8)
#p <- rep(0.3,k)

comb.est <- function(n,k){
  #n=50; k=10
  n.prop <- matrix(0,nrow=k,ncol=2)
  n.prop[,1] <- p
  n.prop[,2] <- 1-p
  tot.seq <- 1:k
  
  ma.ipd.seq <- function(ma.seq){
    #ma.seq <- 1:k
    ipd.seq <- subset(tot.seq,!(tot.seq)%in% ma.seq)
    k1 <- length(ipd.seq)
    k2 <- length(ma.seq)
    n.prop.meta <- n.prop[ma.seq,]
    n.prop.ipd <- n.prop[ipd.seq,]
    var.meta <- n*sum(n.prop.meta[,1]*n.prop.meta[,2])
    b <- sigma.alpha/sigma
    a <- 1+n*b
    var.ipd <- ((n^2)*b/a)*sum(n.prop.ipd[,1]*n.prop.ipd[,2])+
      (n/(a*k1))*(sum(n.prop.ipd[,1])*sum(n.prop.ipd[,2]))
    if(k1==0){
      return(sigma/var.meta)}
    else if (k2==0){
      return(sigma/var.ipd)}
    else{
      return(sigma/(var.meta+var.ipd))}
  }
  #ma.ipd.seq(ma.seq11)
  output <- mapply(ma.ipd.seq,ma.seq.list)
  return(output)
}

k2=0
ma.seq.list <- as.list(data.frame(combn(k,k2)))
no.comb <- ncol(combn(k,k2))
mse_allipd <- comb.est(n,k)

rel.eff <- function(k2){
  ma.seq.list <- as.list(data.frame(combn(k,k2)))
  no.comb <- ncol(combn(k,k2))
  #ma.seq.list <- c(ma.seq.list,list(all.meta=tot.seq),list(all.ipd=NULL))
  mse <- comb.est(n,k)
  rel.eff <- 0.10307346/mse
  # plot(1:no.comb,mse)
  # plot(1:no.comb,rel.eff)
  m <- data.frame(as.character(ma.seq.list),round(mse,4),round(rel.eff,4))
  colnames(m) <- c("ma_seq","mse","rel-eff")
  return(m)
}


k2 <- list(0,2,4)


## all MA mse:0.1295337 rel.eff: 0.7957271

# order <- order(m[,2])
# m_ordered <- cbind(order,m[order,])
# m_ordered <- noquote(m_ordered)
# colnames(m_ordered) <- c("order_no","ma_seq","mse","rel-eff")
# plot(1:no.comb,m_ordered[,3])
# plot(1:no.comb,m_ordered[,4])
# p
# m_ordered
#setwd("C:/Users/Neha Agarwal/Desktop/Neha/Research/meta-analysis/week8/choice2")
#write.csv(m_ordered, file = "result.csv")

dim(m6)
m2 <- m
m4 <- m
m6 <- m
m8 <- m
m2 <- data.frame(m2,comb=1:45)
m4 <- data.frame(m4,comb=1:210)
m6 <- data.frame(m6,comb=1:210)
m8 <- data.frame(m8,comb=1:45)
temp1 <- bind_rows(m8, m6, m4, m2, .id = "id")
colnames(temp1)[colnames(temp1)=="rel-eff"] <- "rel_eff"
temp1$id <- as.factor(temp1$id)
levels(temp1$id) <- c("20% IPD","40% IPD","60% IPD","80% IPD")
  
library(ggplot2)
ggplot(temp1, aes(x=comb,y=rel.eff)) + geom_point() + 
  facet_wrap(~id, scales = "free_x") + 
  geom_hline(yintercept=1,linetype = "dashed") +
  geom_hline(yintercept=0.7957271,color="blue",linetype = "dashed") +
  geom_hline(yintercept=0.9,color="red",linetype = "dashed") +
  labs(x="Different combinations of selected IPD studies",
       y="Relative Efficiency")+ theme_light(base_size = 18)+
  theme(strip.background =element_rect(fill="royalblue1"))+
  theme(strip.text = element_text(colour = 'white'))


m2.maxmin <- rbind(m2[which.max(m2[,3]),2:3],m2[which.min(m2[,3]),2:3])
m4.maxmin <- rbind(m4[which.max(m4[,3]),2:3],m4[which.min(m4[,3]),2:3])
m6.maxmin <- rbind(m6[which.max(m6[,3]),2:3],m6[which.min(m6[,3]),2:3])
m8.maxmin <- rbind(m8[which.max(m8[,3]),2:3],m8[which.min(m8[,3]),2:3])

m2.maxmin <- cbind(m2.maxmin,ma_seq = c("5:6","9:10"))
m4.maxmin <- cbind(m4.maxmin,ma_seq = c("4:7","7:10"))
m6.maxmin <- cbind(m6.maxmin,ma_seq = c("3:8","1,6:10"))
m8.maxmin <- cbind(m8.maxmin,ma_seq = c("2:9","1:2,5:10"))
temp2 <- bind_rows(m8.maxmin,
                   m6.maxmin,
                   m4.maxmin,
                   m2.maxmin,.id = "percent_IPD")
temp2$percent_IPD <- as.factor(temp2$percent_IPD)
levels(temp2$percent_IPD) <- c("20% IPD","40% IPD","60% IPD","80% IPD")


x <- rep(c(1,2,3,4),each=2)
ggplot(temp2, aes(x=x,y=rel.eff,color=percent_IPD)) + geom_point(size=3) + geom_line(lwd=1)+ 
  geom_hline(yintercept=1,linetype = "dashed") +theme_bw() + 
  geom_text(aes(label=ma_seq),hjust=0.40, vjust=0,size=5) +
  geom_hline(yintercept=0.4393,color="blue",linetype = "dashed") +
  labs(x="Maximum and Minimum for different percentage of IPD studies",
       y="Relative Efficiency",size=5)+ theme_light(base_size = 18)+
  theme(strip.background =element_rect(fill="royalblue1"))+
  theme(strip.text = element_text(colour = 'white'))+ 
  scale_color_manual(values=c("blue","red", "skyblue","black"))+ 
  theme(legend.position = c(0.85, 0.2))

order <- order(m8[,2])
m_ordered <- cbind(order,m8[order,])
m_ordered
