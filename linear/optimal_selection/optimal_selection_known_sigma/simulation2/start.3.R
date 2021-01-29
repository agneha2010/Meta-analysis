
optimal.seq.start.3 <- function(k1,prop.trt){
  #k1=4
  #prop.trt=prop.trt.sim[1,]
  ipd.seq.list = as.list(data.frame(t(m3)))
  var = mapply(ma.ipd.seq.expr,ipd.seq.list,MoreArgs = list(prop.trt=prop.trt,k1=3))
  iset = sort(m3[which.min(var),])
  sset = setdiff(seq(1,k), iset)
  sw = sum(y[iset])
  m = sum(prop.trt[iset]*y[iset])/sw
  #v = objfn(prop.trt[iset], y[iset])
  
  for(l in 1:(k1-3)){
    if ((k1-3) <= 0) break
    ind = sset[which.max(y[sset]*(prop.trt[sset]- m)^2/(sw + y[sset]))]
    iset<- union(iset,ind)
    sset = setdiff(sset,ind)
    m = (sw*m + y[ind]*prop.trt[ind])/(sw + y[ind])
    sw = sw + y[ind]
  }
   mat = cbind(sort(iset),ma.ipd.seq.expr(iset,prop.trt,k1))
  return(mat)
}

comb.est.start.3 <- function(prop.trt.sim){
  store = rbind(optimal.seq.start.3(2,prop.trt = prop.trt.sim),
                optimal.seq.start.3(3,prop.trt = prop.trt.sim),
                optimal.seq.start.3(4,prop.trt = prop.trt.sim),
                optimal.seq.start.3(5,prop.trt = prop.trt.sim),
                optimal.seq.start.3(6,prop.trt = prop.trt.sim),
                optimal.seq.start.3(7,prop.trt = prop.trt.sim),
                optimal.seq.start.3(8,prop.trt = prop.trt.sim),
                optimal.seq.start.3(9,prop.trt = prop.trt.sim),
                optimal.seq.start.3(10,prop.trt = prop.trt.sim),
                ma.ipd.seq.expr(NULL,prop.trt.sim,0),
                ma.ipd.seq.expr(1:k,prop.trt.sim,k))
  return(store)
}
result.start.3 <- list()
for(i in 1:nsim){
  print(paste0("simulation ", i))
  result.start.3[[i]] <- round(comb.est.start.3(prop.trt.sim[i,]),4)
}

ipd2 = matrix(NA, nrow=nsim, ncol=7)
ipd3 = matrix(NA, nrow=nsim, ncol=13)
ipd4 = matrix(NA, nrow=nsim, ncol=15)
ipd5 = matrix(NA, nrow=nsim, ncol=17)
ipd6 = matrix(NA, nrow=nsim, ncol=19)
ipd7 = matrix(NA, nrow=nsim, ncol=21)
ipd8 = matrix(NA, nrow=nsim, ncol=23)
ipd9 = matrix(NA, nrow=nsim, ncol=25)
ipd10 = matrix(NA, nrow=nsim, ncol=27)
ipd0 = matrix(NA, nrow=nsim, ncol=2)
ipdk = matrix(NA, nrow=nsim, ncol=2)

check.matches = function(oset1,iset1){
  #oset1 = result[[i]][1:2,1]
  #iset1 = result[[i]][1:2,2]
  k1 = length(oset1)
  t = c()
  for(j in 1:(k1+1)){
    t[j]=ifelse(sum(oset1 %in% iset1)==(j-1),1,0)
  }
  return(t)
}

for(i in 1:nsim){
  ipd2[i,] <- c(check.matches(result[[i]][1:2,1],result[[i]][1:2,2]),result[[i]][1,3:6])
  ipd3[i,] <- c(check.matches(result[[i]][3:5,1],result[[i]][3:5,2]),result[[i]][3,3:6],
                check.matches(result[[i]][3:5,1],result.start.3[[i]][4:6,1]),result.start.3[[i]][4,2])
  ipd4[i,] <- c(check.matches(result[[i]][6:9,1],result[[i]][6:9,2]),result[[i]][6,3:6],
                check.matches(result[[i]][6:9,1],result.start.3[[i]][7:10,1]),result.start.3[[i]][7,2])
  ipd5[i,] <- c(check.matches(result[[i]][10:14,1],result[[i]][10:14,2]),result[[i]][10,3:6],
                check.matches(result[[i]][10:14,1],result.start.3[[i]][11:15,1]),result.start.3[[i]][11,2])
  ipd6[i,] <- c(check.matches(result[[i]][15:20,1],result[[i]][15:20,2]),result[[i]][15,3:6],
                check.matches(result[[i]][15:20,1],result.start.3[[i]][16:21,1]),result.start.3[[i]][16,2])
  ipd7[i,] <- c(check.matches(result[[i]][21:27,1],result[[i]][21:27,2]),result[[i]][21,3:6],
                check.matches(result[[i]][21:27,1],result.start.3[[i]][22:28,1]),result.start.3[[i]][22,2])
  ipd8[i,] <- c(check.matches(result[[i]][28:35,1],result[[i]][28:35,2]),result[[i]][28,3:6],
                check.matches(result[[i]][28:35,1],result.start.3[[i]][29:36,1]),result.start.3[[i]][29,2])
  ipd9[i,] <- c(check.matches(result[[i]][36:44,1],result[[i]][36:44,2]),result[[i]][38,3:6],
                check.matches(result[[i]][36:44,1],result.start.3[[i]][37:45,1]),result.start.3[[i]][37,2])
  ipd10[i,] <- c(check.matches(result[[i]][45:54,1],result[[i]][45:54,2]),result[[i]][45,3:6],
                 check.matches(result[[i]][45:54,1],result.start.3[[i]][46:55,1]),result.start.3[[i]][46,2])
  ipd0[i,] <- c(result[[i]][55,1],result.start.3[[i]][56,1])
  ipdk[i,] <- c(result[[i]][56,1],result.start.3[[i]][57,1])
}

### number of places in which they match == 0,1,...,10
colnames(ipd2) = c(paste0("match_count_",0:2),"min_var_oset","var_iset",
                   "max_var_oset","ratio_obj_iset_oset")
colnames(ipd3) = c(paste0("match_count_",0:3),"min_var_oset","var_iset",
                   "max_var_oset","ratio_obj_iset_oset",paste0("match_count_start_3",0:3),"var_iset_start_3")
colnames(ipd4) = c(paste0("match_count_",0:4),"min_var_oset","var_iset",
                   "max_var_oset","ratio_obj_iset_oset",paste0("match_count_start_3",0:4),"var_iset_start_3")
colnames(ipd5) = c(paste0("match_count_",0:5),"min_var_oset","var_iset",
                   "max_var_oset","ratio_obj_iset_oset",paste0("match_count_start_3",0:5),"var_iset_start_3")
colnames(ipd6) = c(paste0("match_count_",0:6),"min_var_oset","var_iset",
                   "max_var_oset","ratio_obj_iset_oset",paste0("match_count_start_3",0:6),"var_iset_start_3")
colnames(ipd7) = c(paste0("match_count_",0:7),"min_var_oset","var_iset",
                   "max_var_oset","ratio_obj_iset_oset",paste0("match_count_start_3",0:7),"var_iset_start_3")
colnames(ipd8) = c(paste0("match_count_",0:8),"min_var_oset","var_iset",
                   "max_var_oset","ratio_obj_iset_oset",paste0("match_count_start_3",0:8),"var_iset_start_3")
colnames(ipd9) = c(paste0("match_count_",0:9),"min_var_oset","var_iset",
                   "max_var_oset","ratio_obj_iset_oset",paste0("match_count_start_3",0:9),"var_iset_start_3")
colnames(ipd10) = c(paste0("match_count_",0:10),"min_var_oset","var_iset",
                    "max_var_oset","ratio_obj_iset_oset",paste0("match_count_start_3",0:10),"var_iset_start_3")
count_matches = cbind(colSums(ipd10[,1:11]),colSums(ipd9[,1:10]),colSums(ipd8[,1:9]),
                      colSums(ipd7[,1:8]),colSums(ipd6[,1:7]),colSums(ipd5[,1:6]),
                      colSums(ipd4[,1:5]),colSums(ipd3[,1:4]),colSums(ipd2[,1:3]))
colnames(count_matches) = paste0("k1_",10:2)

count_matches_start_3 = cbind(colSums(ipd10[,16:26]),colSums(ipd9[,15:24]),colSums(ipd8[,14:22]),
                      colSums(ipd7[,13:20]),colSums(ipd6[,12:18]),colSums(ipd5[,11:16]),
                      colSums(ipd4[,10:14]),colSums(ipd3[,9:12]))
colnames(count_matches_start_3) = paste0("k1_",10:3)


mean_ratio = c(mean(ipd10[,15]),mean(ipd9[,14]),mean(ipd8[,13]),
               mean(ipd7[,12]),mean(ipd6[,11]),mean(ipd5[,10]),
               mean(ipd4[,9]),mean(ipd3[,8]),mean(ipd2[,7]))
count_matches = rbind(count_matches,mean_ratio)

mean_ratio_var_iset_oset = c(mean(ipd10[,13]/ipd10[,12]),mean(ipd9[,12]/ipd9[,11]),
                         mean(ipd8[,11]/ipd8[,10]),mean(ipd7[,10]/ipd7[,9]),
                         mean(ipd6[,9]/ipd6[,8]),mean(ipd5[,8]/ipd5[,7]),
                         mean(ipd4[,7]/ipd4[,6]),mean(ipd3[,6]/ipd3[,5]),
                         mean(ipd2[,5]/ipd2[,4]))
count_matches = rbind(count_matches,mean_ratio_var_iset_oset)

mean_ratio_var_iset_start_3_oset = c(mean(ipd10[,27]/ipd10[,12]),mean(ipd9[,25]/ipd9[,11]),
                         mean(ipd8[,23]/ipd8[,10]),mean(ipd7[,21]/ipd7[,9]),
                         mean(ipd6[,19]/ipd6[,8]),mean(ipd5[,17]/ipd5[,7]),
                         mean(ipd4[,15]/ipd4[,6]),mean(ipd3[,13]/ipd3[,5]))
count_matches_start_3 = rbind(count_matches_start_3,mean_ratio_var_iset_start_3_oset)

# write.csv(count_matches_start_3,"count_matches_start_3.csv")
# write.csv(count_matches,"count_matches.csv")
# write.csv(ipd2,"ipd2.csv")
# write.csv(ipd3,"ipd3.csv")
# write.csv(ipd4,"ipd4.csv")
# write.csv(ipd5,"ipd5.csv")
# write.csv(ipd6,"ipd6.csv")
# write.csv(ipd7,"ipd7.csv")
# write.csv(ipd8,"ipd8.csv")
# write.csv(ipd9,"ipd9.csv")
# write.csv(ipd10,"ipd10.csv")
# write.csv(ipd0,"ipd0.csv")
# write.csv(ipdk,"ipdk.csv")
