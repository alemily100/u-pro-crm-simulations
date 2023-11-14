setwd("M:/PhD/Trial Designs/ProTox")
source("R/functions.R")
library(dfcrm)
########################################################################################################
no.sim<-5000


###########################Run block before every simulation###############################################################################

u<- c(0.06, 0.14, 0.25, 0.38, 0.50)
v<- c(0.10, 0.21, 0.35, 0.49, 0.61)
target_c<- 0.25
target_p<- 0.35


#Simulation scenarios 
true_tox_clin<- rbind.data.frame(c(0.05, 0.05, 0.25, 0.4, 0.55), c(0.05, 0.25, 0.40, 0.55, 0.7), c(0.01, 0.02, 0.05, 0.10, 0.25),
                                 c(0.02, 0.05, 0.1, 0.25, 0.4), c(0.05, 0.1, 0.16, 0.25, 0.4), c(0.05, 0.18, 0.2, 0.25, 0.4), 
                                 c(0.01, 0.05, 0.1, 0.16, 0.25),c(0.25, 0.40, 0.55, 0.7, 0.8), c(0.45, 0.50, 0.55, 0.7, 0.8))
dimnames(true_tox_clin)[[2]]<- c( "1","2", "3", "4", "5")
dimnames(true_tox_clin)[[1]]<- c("Sc 1", "Sc 2", "Sc 3", "Sc 4", "Sc 5", "Sc 6", "Sc 7", "Sc 8", "Sc 9")                        

true_tox_pro<- rbind.data.frame(c(0.17, 0.18, 0.35, 0.50, 0.65), c(0.1, 0.15, 0.35, 0.5, 0.65), c(0.04, 0.09, 0.17, 0.2, 0.35),
                                c(0.09, 0.17, 0.2, 0.35, 0.5), c(0.05, 0.2, 0.35, 0.5, 0.65), c(0.17, 0.35, 0.5, 0.65, 0.8),
                                c(0.04, 0.05, 0.2, 0.35, 0.5),c(0.35, 0.5, 0.65, 0.8, 0.85), c(0.55, 0.6, 0.65, 0.8, 0.85))
dimnames(true_tox_pro)[[2]]<- c( "1","2", "3", "4", "5")
dimnames(true_tox_pro)[[1]]<- c("Sc 1", "Sc 2", "Sc 3", "Sc 4", "Sc 5", "Sc 6", "Sc 7","Sc 8", "Sc 9")                        

no_enrolled<-3
phi<-0.1
alpha<- 15
no.dosages<- 5
sample<- 39
stop_prob<-0.8
target_stop<- 0.5
all_sc<-NULL
a<- 0.35
b<- 0.65
############################original PRO-CRM design######################################################################################
############################phi<- 0.1####################################################################
prop_mtd<-c()
prop_overdosed<-c()
#start<- Sys.time()
for (i in 1:9){
  cl <- makeCluster(2)
  clusterSetRNGStream(cl, 2115)
  invisible(clusterEvalQ(cl,{
    setwd("M:/PhD/Trial Designs/ProTox")
    source("R/functions.R")
  }))
  
  clin_true_tox<- as.numeric(true_tox_clin[i,])
  pat_true_tox<- as.numeric(true_tox_pro[i,])
  
  clusterExport(cl, c("u", "v", "no.sim", "clin_true_tox","sample","all_sc","target_stop",
                      "pat_true_tox", "target_c", "target_p", "no_enrolled", "phi", "no.dosages", "stop_prob", "prop_mtd", "prop_overdosed", "a", "b"))
  
  collect<-parSapply(cl,rep(sample, times=no.sim), function(k) trial_sim_original(u,v,k, no_enrolled, phi, clin_true_tox, pat_true_tox, no.dosages,target_c, target_p, stop_prob, a, b, target_stop)) 
  dlt<- collect[1,]
  prop_mtd[i]<- mean(collect[2,], na.rm=TRUE)
  prop_overdosed[i]<- mean(collect[3,], na.rm=TRUE)
  stopCluster(cl)
  prob_dlt<- c(sum(dlt==1,na.rm=TRUE)/no.sim, sum(dlt==2,na.rm=TRUE)/no.sim, sum(dlt==3,na.rm=TRUE)/no.sim, sum(dlt==4,na.rm=TRUE)/no.sim, sum(dlt==5,na.rm=TRUE)/no.sim, sum(is.na(dlt))/no.sim)
  sc<-rbind.data.frame(c(clin_true_tox, NA), c(pat_true_tox, NA), prob_dlt)
  dimnames(sc)[[2]]<- c( "1","2", "3", "4", "5", "stop")
  dimnames(sc)[[1]]<- c("Probability of clinician DLT", "Probability of patient DLT", "Probability recommended as MTD")
  all_sc<-rbind.data.frame(all_sc, sc)
  write.csv(all_sc, paste0(alpha, "_original_sim_",phi,".csv"))
}

#time<- Sys.time()-start
#print(time)
allocated_mtd<-rbind.data.frame(prop_mtd)
dimnames(allocated_mtd)[[2]]<-paste0("scenario",c(1:9))
write.csv(allocated_mtd, paste0(alpha, "_original_allocated_",phi,".csv"))

overdosed<-rbind.data.frame(prop_overdosed)
dimnames(overdosed)[[2]]<-paste0("scenario",c(1:9))
write.csv(overdosed, paste0(alpha, "_original_overdosed_",phi,".csv"))

############################utility PRO-CRM design######################################################################################
################################phi = 0.1 ####################################################
prop_mtd<-c()
prop_overdosed<-c()
all_sc<-NULL
start<- Sys.time()
for (i in 1:9){
  cl <- makeCluster(2)
  clusterSetRNGStream(cl, 2115)
  invisible(clusterEvalQ(cl,{
    setwd("M:/PhD/Trial Designs/ProTox")
    source("R/functions.R")
  }))
  
  clin_true_tox<- as.numeric(true_tox_clin[i,])
  pat_true_tox<- as.numeric(true_tox_pro[i,])
  
  clusterExport(cl, c("u", "v", "alpha", "no.sim", "clin_true_tox","sample","all_sc","target_stop",
                      "pat_true_tox", "target_c", "target_p", "no_enrolled", "phi", "no.dosages", "stop_prob","prop_mtd", "prop_overdosed","a", "b"))
  collect<-parSapply(cl,rep(sample, times=no.sim), function(k) trial_sim_utility(u,v,alpha,k,no_enrolled, phi, clin_true_tox, pat_true_tox, no.dosages,target_c, target_p, stop_prob, a, b, target_stop)) 
  dlt<- collect[1,]
  prop_mtd[i]<- mean(collect[2,], na.rm=TRUE)
  prop_overdosed[i]<- mean(collect[3,], na.rm=TRUE)
  stopCluster(cl)
  prob_dlt<- c(sum(dlt==1,na.rm=TRUE)/no.sim, sum(dlt==2,na.rm=TRUE)/no.sim, sum(dlt==3,na.rm=TRUE)/no.sim, sum(dlt==4,na.rm=TRUE)/no.sim, sum(dlt==5,na.rm=TRUE)/no.sim, sum(is.na(dlt))/no.sim)
  sc<-rbind.data.frame(c(clin_true_tox, NA), c(pat_true_tox, NA), prob_dlt)
  dimnames(sc)[[2]]<- c( "1","2", "3", "4", "5", "stop")
  dimnames(sc)[[1]]<- c("Probability of clinician DLT", "Probability of patient DLT", "Probability recommended as MTD")
  all_sc<-rbind.data.frame(all_sc, sc)
  write.csv(all_sc, paste0(alpha, "_utility_sim_",phi,".csv"))
}

allocated_mtd<-rbind.data.frame(prop_mtd)
dimnames(allocated_mtd)[[2]]<-paste0("scenario",c(1:9))
write.csv(allocated_mtd, paste0(alpha, "_utility_allocated_",phi,".csv"))

overdosed<-rbind.data.frame(prop_overdosed)
dimnames(overdosed)[[2]]<-paste0("scenario",c(1:9))
write.csv(overdosed, paste0(alpha, "_utility_overdosed_",phi,".csv"))

############################original PRO-CRM design######################################################################################
############################phi<- 0.9####################################################################
phi<- 0.9
all_sc<- NULL
prop_mtd<-c()
prop_overdosed<-c()
#start<- Sys.time()
for (i in 1:9){
  cl <- makeCluster(2)
  clusterSetRNGStream(cl, 2115)
  invisible(clusterEvalQ(cl,{
    setwd("M:/PhD/Trial Designs/ProTox")
    source("R/functions.R")
  }))
  
  clin_true_tox<- as.numeric(true_tox_clin[i,])
  pat_true_tox<- as.numeric(true_tox_pro[i,])
  
  clusterExport(cl, c("u", "v", "no.sim", "clin_true_tox","sample","all_sc","target_stop",
                      "pat_true_tox", "target_c", "target_p", "no_enrolled", "phi", "no.dosages", "stop_prob", "prop_mtd", "prop_overdosed","a", "b"))
  
  collect<-parSapply(cl,rep(sample, times=no.sim), function(k) trial_sim_original(u,v,k, no_enrolled, phi, clin_true_tox, pat_true_tox, no.dosages,target_c, target_p, stop_prob,a,b, target_stop)) 
  dlt<- collect[1,]
  prop_mtd[i]<- mean(collect[2,], na.rm=TRUE)
  prop_overdosed[i]<- mean(collect[3,], na.rm=TRUE)
  stopCluster(cl)
  prob_dlt<- c(sum(dlt==1,na.rm=TRUE)/no.sim, sum(dlt==2,na.rm=TRUE)/no.sim, sum(dlt==3,na.rm=TRUE)/no.sim, sum(dlt==4,na.rm=TRUE)/no.sim, sum(dlt==5,na.rm=TRUE)/no.sim, sum(is.na(dlt))/no.sim)
  sc<-rbind.data.frame(c(clin_true_tox, NA), c(pat_true_tox, NA), prob_dlt)
  dimnames(sc)[[2]]<- c( "1","2", "3", "4", "5", "stop")
  dimnames(sc)[[1]]<- c("Probability of clinician DLT", "Probability of patient DLT", "Probability recommended as MTD")
  all_sc<-rbind.data.frame(all_sc, sc)
  write.csv(all_sc, paste0(alpha, "_original_sim_",phi,".csv"))
}

#time<- Sys.time()-start
#print(time)
allocated_mtd<-rbind.data.frame(prop_mtd)
dimnames(allocated_mtd)[[2]]<-paste0("scenario",c(1:9))
write.csv(allocated_mtd, paste0(alpha, "_original_allocated_",phi,".csv"))

overdosed<-rbind.data.frame(prop_overdosed)
dimnames(overdosed)[[2]]<-paste0("scenario",c(1:9))
write.csv(overdosed, paste0(alpha, "_original_overdosed_",phi,".csv"))



############################utility PRO-CRM design######################################################################################
################################phi = 0.1 ####################################################
prop_mtd<-c()
prop_overdosed<-c()
all_sc<-NULL
start<- Sys.time()
for (i in 1:9){
  cl <- makeCluster(2)
  clusterSetRNGStream(cl, 2115)
  invisible(clusterEvalQ(cl,{
    setwd("M:/PhD/Trial Designs/ProTox")
    source("R/functions.R")
  }))
  
  clin_true_tox<- as.numeric(true_tox_clin[i,])
  pat_true_tox<- as.numeric(true_tox_pro[i,])
  
  clusterExport(cl, c("u", "v", "alpha", "no.sim", "clin_true_tox","sample","all_sc","target_stop",
                      "pat_true_tox", "target_c", "target_p", "no_enrolled", "phi", "no.dosages", "stop_prob","prop_mtd", "prop_overdosed","a", "b"))
  collect<-parSapply(cl,rep(sample, times=no.sim), function(k) trial_sim_utility(u,v,alpha,k,no_enrolled, phi, clin_true_tox, pat_true_tox, no.dosages,target_c, target_p, stop_prob, a,b, target_stop)) 
  dlt<- collect[1,]
  prop_mtd[i]<- mean(collect[2,], na.rm=TRUE)
  prop_overdosed[i]<- mean(collect[3,], na.rm=TRUE)
  stopCluster(cl)
  prob_dlt<- c(sum(dlt==1,na.rm=TRUE)/no.sim, sum(dlt==2,na.rm=TRUE)/no.sim, sum(dlt==3,na.rm=TRUE)/no.sim, sum(dlt==4,na.rm=TRUE)/no.sim, sum(dlt==5,na.rm=TRUE)/no.sim, sum(is.na(dlt))/no.sim)
  sc<-rbind.data.frame(c(clin_true_tox, NA), c(pat_true_tox, NA), prob_dlt)
  dimnames(sc)[[2]]<- c( "1","2", "3", "4", "5", "stop")
  dimnames(sc)[[1]]<- c("Probability of clinician DLT", "Probability of patient DLT", "Probability recommended as MTD")
  all_sc<-rbind.data.frame(all_sc, sc)
  write.csv(all_sc, paste0(alpha,"_utility_sim_",phi,".csv"))
}

allocated_mtd<-rbind.data.frame(prop_mtd)
dimnames(allocated_mtd)[[2]]<-paste0("scenario",c(1:9))
write.csv(allocated_mtd, paste0(alpha,"_utility_allocated_",phi,".csv"))

overdosed<-rbind.data.frame(prop_overdosed)
dimnames(overdosed)[[2]]<-paste0("scenario",c(1:9))
write.csv(overdosed, paste0(alpha, "_utility_overdosed_",phi,".csv"))
