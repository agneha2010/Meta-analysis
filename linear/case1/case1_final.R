
library(mvtnorm)
library(dplyr)
library(Matrix)

k <- 10
n <- 10
beta <- 1.5 ;alpha <- 0.5; sigma <- 5
sigma.alpha <- 0.025
true.theta <- c(alpha,beta)
set.seed(100)
alpha.k <- rnorm(k,alpha,sqrt(sigma.alpha))

p <- c(0.1, 0.2, 0.3, 0.3, 0.3, 0.5, 0.6, 0.6, 0.8, 0.8)
n.prop <- matrix(0,nrow=k,ncol=2)
n.prop[,1] <- p
n.prop[,2] <- 1-p
tot.seq <- 1:k
  
ma.ipd.seq <- function(ma.seq){

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
output <- mapply(ma.ipd.seq,ma.seq.list)

k2=0
ma.seq.list <- as.list(data.frame(combn(k,k2)))
mse_allipd <- ma.ipd.seq(NULL)
mse_allmeta <- ma.ipd.seq(1:k)

(mse_allmeta-mse_allipd)/mse_allipd

rel.eff <- function(k2.par){
  #k2.par = 4
  ma.seq.list <- as.list(data.frame(combn(k,k2.par)))
  no.comb <- ncol(combn(k,k2.par))
  mse <- mapply(ma.ipd.seq,ma.seq.list)
  rel.eff <- mse_allipd/mse
  m <- data.frame(as.character(ma.seq.list),round(mse,4),round(rel.eff,4))
  colnames(m) <- c("ma_seq","mse","rel-eff")
  return(m)
}

dim(m6)
m2 <- rel.eff(2) #80% IPD
m4 <- rel.eff(4)
m6 <- rel.eff(6)
m8 <- rel.eff(8) #20% IPD
m2 <- data.frame(m2,comb=1:nrow(m2))
m4 <- data.frame(m4,comb=1:nrow(m4))
m6 <- data.frame(m6,comb=1:nrow(m6))
m8 <- data.frame(m8,comb=1:nrow(m8))
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
