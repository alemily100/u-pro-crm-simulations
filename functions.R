###########################load_packages#######################
library(Rlab)
library(MASS)
library(parallel)


###############################like_c#############################
#Description: evaluate log-likelihood for C-DLT endpoint

#Input:
#gamma - estimate of parameter 
#data_frame - data frame with column 1 containing the patient ID, column 2 containing therr dosage and column 3 containing their C-DLT outcome 
#u - skeleton for C-DLT rate 

#Output: Numerical; log likelihood for parameter using C-DLT endpoint 

like_c<- function(gamma, data_frame,u){
  l<-0
  outcome<-data_frame[,3]
  dose<- data_frame[,2]
  for(i in 1:nrow(data_frame)){
    l<-l+(as.numeric(outcome[i]==0))*log(1-u[dose[i]]^gamma)+(as.numeric(outcome[i]==1)*gamma)*log(u[dose[i]])
  }
  l
}

###############################like_p#############################
#Description: evaluate log-likelihood for C-DLT endpoint

#Input:
#eta - estimate of parameter 
#data_frame - data frame with column 1 containing the patient ID, column 2 containing their dosage and column 3 containing their P-DLT outcome 
#v - skeleton for P-DLT rate 

#Output: Numerical; log likelihood for parameter using P-DLT endpoint  

like_p<- function(eta, data_frame,v){
  l<-0
  outcome<-data_frame[,4]
  dose<- data_frame[,2]
  for(i in 1:nrow(data_frame)){
    l<-l+(as.numeric(outcome[i]==0))*log(1-v[dose[i]]^eta)+(as.numeric(outcome[i]==1)*eta)*log(v[dose[i]])
  }
  l
}

###########################time_to_dlt#######################################
#Description: Simulates correlated C-DLT and P-DLT as per the Clayton model

#Input: 
#dose - dose for patient i
#clin_tox - vector of true probability of C-DLT for each of the 5 dosages 
#pat_tox - vector of true probability of P-DLT for each of the 5 dosages 
#phi - correlation parameter

#Output: Vector; time to C-DLT and P-DLT for patient i 

time_to_dlt<- function(dose,clin_tox, pat_tox, phi){
  lambda_c<- -log(1-clin_tox[dose])
  lambda_p<- -log(1-pat_tox[dose])
  u1<-runif(1)
  time_c<- (-log(u1))/lambda_c
  u2<-runif(1)
  a<- u2/(u1^(-(phi+1)/phi))
  b<- u1^(-1/phi)-1
  time_p<- (phi/lambda_p)*log(a^(1/(-phi-1))-b)
  times<- c(time_c, time_p)
  return(times)
}

#########################truncated_matrix##################################
#Description: Function to generate data frame of all patients, their dosages and DLT endpoints

#Input: 
#no_patients - the number of new patients enrolled in the trial at cycle t 
#matrix_time_to_toxicity - matrix of the new patient's C-DLT and P-DLT times, as generated using `time_to_dlt`
#previous_matrix - previous matrix which contains all information up until just before t
#dose - dose given to the new patient enrolled in the trial at time t

#Output: Data frame; containing all patients, their dosages and DLT endpoints at cycle t

truncated_matrix<- function(no_patients, matrix_time_to_toxicity, previous_matrix, dose){
  new_mat<- matrix(nrow=no_patients, ncol=4)
  new_mat[,1]<-(max(previous_matrix[,1])+1):(max(previous_matrix[,1])+no_patients)
  new_mat[,2]<-rep(dose, times=no_patients)
  new_mat[,3]<-as.numeric(matrix_time_to_toxicity[,1]<=1)
  new_mat[,4]<-as.numeric(matrix_time_to_toxicity[,2]<=1)
  return(rbind(previous_matrix, new_mat))
}

####################stopping_rule###################################
#Description: Assess whether a trial should be stopped for safety pre-maturely 

#Input: 
#matrix: data frame containing all patients, their dosages and DLT endpoints at cycle t
#a - shape parameter 
#b - shape parameter 
#target - \tilde{\pi}_{CP}
#min_p - \tau

#Output: Text; confirming whether a trial should be stopped prematurely ("stop") or not ("continue") 

stopping_rule<-function(matrix, a,b, target, min_p){
  m_1<- matrix[matrix[,2]==1,]
  sum_dlt<- sum((m_1[,3] | m_1[,4]))
  a_tilde<- a + sum_dlt
  b_tilde<- b + nrow(m_1) - sum_dlt
  ifelse(pbeta(target, a_tilde, b_tilde, lower.tail = FALSE)>min_p,
         "stop", "continue")
}

#################################f##############################################
#Description: utility curve function

#Input:
#a - alpha 
#pat_prob - Estimated probability of P-DLT 
#clin_prob - Estimated probability of C-DLT
#pat_target - target P-DLT rate 
#clin_target - target C-DLT rate

#Output: Numerical value; calculating LHS of Equation 1 in the manuscript

f<-function(a,pat_prob,clin_prob, pat_target, clin_target){
  1-((pat_prob/pat_target)^a + (clin_prob/clin_target)^a)^(1/a)
}

###############################distance########################################
#Description - Euclidean distance from estimated C-DLT and P-DLT for a dose and a specified poinr on the utilty curve 

#Input: 
#alpha - alpha of utility 
# x - point on the utility curve used to calculate distance 
#pred_x - estimated P-DLT rate for dose
#pred_y - estimated C-DLT rate for dose 
# patient_target_DLT - target P-DLT rate
# clinician_target_DLT - target C-DLT rate

#Output: Numerical; The Euclidean distance between estimated (P-DLT rate, C-DLT rate) and (x,f(x)) where f defines the utility curve 

distance<-function(alpha, x, pred_x, pred_y,patient_target_DLT,clinician_target_DLT){
  f_x<- uniroot(f, a=alpha, pat_prob=x, pat_target=patient_target_DLT, clin_target=clinician_target_DLT,interval= c(1.e-14, 1e04),
                extendInt="yes")$root
  d<-sqrt((x-pred_x)^2+(f_x-pred_y)^2)
  return(d)
}

##############################trial_sim_original####################################
#Description: Generate MTD selection for one generated trial under original PRO-CRM decision criteria 

#u_skeleton - skeleton for C-DLT
#v_skeleton - skeleton for P-DLT
#sample - overall sample size of trial 
#no_enrolled - number of patients enrolled after each interim analysis 
#phi - correlation parameter for DLT time generation 
#true_tox_c - true C-DLT rate under the simulation scenario 
#true_tox_p - true P-DLT rate under simulation scenario 
#number_dosages - number of dosages under investigation 
#target_clin - target C-DLT rate
#target_pat - target P-DLT rate 
#stopping_rule_prob<- \tau
#a - shape parameter for beta distribution 
#b - shape parameter for beta distribution 
#target<- \tilde{\pi}_{CP}

#Ouput: Vector; containing recommended MTD, proportion of patients assigned the MTD, and the proporion of patients overdosed. 

trial_sim_original<- function(u_skeleton, v_skeleton, sample,no_enrolled, phi, true_tox_c, true_tox_p, number_dosages, target_clin, target_pat,stopping_rule_prob, a, b, target){
  current_dose<-1
  val<-t(sapply(rep(current_dose, times=3), function (k) time_to_dlt(k,true_tox_c, true_tox_p, phi)))
  M<- matrix(nrow=3, ncol=4)
  M[,1]<-1:3
  M[,2]<-rep(current_dose, times=3)
  M[,3]<-val[,1]<=1
  M[,4]<-val[,2]<=1
  while(sum(M[,3], M[,4])==0 && nrow(M)<sample){
    #use rule based design when there is no heterogeneity in either endpoint 
    if(stopping_rule(M, a,b, target, stopping_rule_prob)=="stop"){
      current_dose<- NA
      break
    }
    ifelse(current_dose<number_dosages,current_dose<-current_dose+1, current_dose<-number_dosages)
    val<-t(sapply(rep(current_dose, times=3), function (k) time_to_dlt(k,true_tox_c, true_tox_p, phi)))
    M<-truncated_matrix(no_enrolled, val, M, current_dose)
  }
  while(xor(sum(M[,3])>0,sum(M[,4])>0) & nrow(M)<sample){
    #use rule based design and likelihood CRM when there is heterogeneity in one endpoint only  
    if(stopping_rule(M, 0.1,0.9, target, stopping_rule_prob)=="stop"){
      current_dose<- NA
      break
    }
    if(sum(M[,3])>0){
      est_c<-optimise(like_c,c(0,10), data_frame=M,u=u_skeleton, maximum = TRUE)
      p_c<-u_skeleton^est_c$maximum 
      current_dose<-min(which.min(abs(p_c-target_clin)), current_dose+1)
      val<-t(sapply(rep(current_dose, times=no_enrolled), function (k) time_to_dlt(k,true_tox_c, true_tox_p, phi)))
      M<-truncated_matrix(no_enrolled, val, M, current_dose)
    }
    else{
      est_p<-optimise(like_p, c(0,10), data_frame=M,v=v_skeleton,maximum = TRUE)
      p_p <- v_skeleton^est_p$maximum
      current_dose<-min(which.min(abs(p_p-target_pat)), current_dose+1)
      val<-t(sapply(rep(current_dose, times=no_enrolled), function (k) time_to_dlt(k,true_tox_c, true_tox_p, phi)))
      M<-truncated_matrix(no_enrolled, val, M, current_dose)
    }
  }
  while(sum(M[,3])>0&sum(M[,4])>0&nrow(M)<sample){
    #use likelihood CRM when there is heterogeneity in both endpoints   
    if(stopping_rule(M, a,b, target, stopping_rule_prob)=="stop"){
      current_dose<- NA
      break
    }
    est_c<-optimise(like_c,c(0,10), u=u_skeleton, data_frame=M, maximum = TRUE)
    est_p<-optimise(like_p,c(0,10), v=v_skeleton, data_frame=M, maximum = TRUE)
    #dose-finding decision 
    p_c<-u_skeleton^est_c$maximum 
    p_p<-v_skeleton^est_p$maximum
    rec_dose<-min(which.min(abs(p_c-target_clin)), which.min(abs(p_p-target_pat)))
    ifelse(rec_dose>max(M[,2]), current_dose<-max(M[,2])+1, current_dose<-rec_dose)
    val<-t(sapply(rep(current_dose, times=no_enrolled), function (k) time_to_dlt(k,true_tox_c, true_tox_p, phi)))
    M<-truncated_matrix(no_enrolled, val, M, current_dose)
  }
   est_c<-optimise(like_c,c(0,10), u=u_skeleton, data_frame=M, maximum = TRUE)
   est_p<-optimise(like_p,c(0,10), v=v_skeleton, data_frame=M, maximum = TRUE)
   #dose-finding decision 
   p_c<-u_skeleton^est_c$maximum 
   p_p<-v_skeleton^est_p$maximum
   rec_dose<-min(which.min(abs(p_c-target_clin)), which.min(abs(p_p-target_pat)))
  ifelse(rec_dose>max(M[,2]), current_dose<-max(M[,2])+1, current_dose<-rec_dose)
  unlink(file.path("C:/Users/ealger/AppData/Local/Temp", "Rtmp*"), recursive = T)
  optimal_dose<- min(which.min(abs(true_tox_c-target_clin)), which.min(abs(true_tox_p-target_pat)))
  mtd_assigned<-sum(M[,2]==optimal_dose)/sample
  overdosed<- sum(M[,2]>optimal_dose)/sample
  return(c(current_dose, mtd_assigned,overdosed))
}



##############################trial_sim_utility####################################
#Description: Generate MTD selection for one generated trial under U-PRO-CRM decision criteria 

#u_skeleton - skeleton for C-DLT
#v_skeleton - skeleton for P-DLT
#alpha - alpha used to define the utility curve  
#sample - overall sample size of trial 
#no_enrolled - number of patients enrolled after each interim analysis 
#phi - correlation parameter for DLT time generation 
#true_tox_c - true C-DLT rate under the simulation scenario 
#true_tox_p - true P-DLT rate under simulation scenario 
#number_dosages - number of dosages under investigation 
#target_clin - target C-DLT rate
#target_pat - target P-DLT rate 
#stopping_rule_prob<- \tau
#a - shape parameter for beta distribution 
#b - shape parameter for beta distribution 
#target<- \tilde{\pi}_{CP}

#Ouput: Vector; containing recommended MTD, proportion of patients assigned the MTD, and the proportion of patients overdosed. 

trial_sim_utility<- function(u_skeleton, v_skeleton, alpha, sample,no_enrolled, phi, true_tox_c, true_tox_p, number_dosages, target_clin, target_pat, stopping_rule_prob, a,b, target){
  current_dose<-1
  val<-t(sapply(rep(current_dose, times=3), function (k) time_to_dlt(k,true_tox_c, true_tox_p, phi)))
  M<- matrix(nrow=3, ncol=4)
  M[,1]<-1:3
  M[,2]<-rep(current_dose, times=3)
  M[,3]<-val[,1]<=1
  M[,4]<-val[,2]<=1
  while(sum(M[,3], M[,4])==0 && nrow(M)<sample){
    #use rule based design when there is no heterogeneity in either endpoint 
    if(stopping_rule(M, a,b, target, stopping_rule_prob)=="stop"){
      current_dose<- NA
      break
    }
    ifelse(current_dose<number_dosages,current_dose<-current_dose+1, current_dose<-number_dosages)
    val<-t(sapply(rep(current_dose, times=3), function (k) time_to_dlt(k,true_tox_c, true_tox_p, phi)))
    M<-truncated_matrix(no_enrolled, val, M, current_dose)
  }
  while(xor(sum(M[,3])>0,sum(M[,4])>0) & nrow(M)<sample){
    #use rule based design and likelihood CRM when there is heterogeneity in one endpoint only  
    if(stopping_rule(M, 0.1,0.9, target, stopping_rule_prob)=="stop"){
      current_dose<- NA
      break
    }
    if(sum(M[,3])>0){
      est_c<-optimise(like_c,c(0,10), data_frame=M,u=u_skeleton, maximum = TRUE)
      p_c<-u_skeleton^est_c$maximum 
      current_dose<-min(which.min(abs(p_c-target_clin)), current_dose+1)
      val<-t(sapply(rep(current_dose, times=no_enrolled), function (k) time_to_dlt(k,true_tox_c, true_tox_p, phi)))
      M<-truncated_matrix(no_enrolled, val, M, current_dose)
    }
    else{
      est_p<-optimise(like_p, c(0,10), data_frame=M,v=v_skeleton,maximum = TRUE)
      p_p <- v_skeleton^est_p$maximum
      current_dose<-min(which.min(abs(p_p-target_pat)), current_dose+1)
      val<-t(sapply(rep(current_dose, times=no_enrolled), function (k) time_to_dlt(k,true_tox_c, true_tox_p, phi)))
      M<-truncated_matrix(no_enrolled, val, M, current_dose)
    }
  }
  while(sum(M[,3])>0&sum(M[,4])>0&nrow(M)<sample){
    #use likelihood CRM when there is heterogeneity in both endpoints   
    if(stopping_rule(M, a,b, target, stopping_rule_prob)=="stop"){
      current_dose<- NA
      break
    }
    est_c<-optimise(like_c,c(0,10), u=u_skeleton, data_frame=M, maximum = TRUE)
    est_p<-optimise(like_p,c(0,10), v=v_skeleton, data_frame=M, maximum = TRUE)
    #dose-finding decision 
    p_c<-u_skeleton^est_c$maximum 
    p_p<-v_skeleton^est_p$maximum
    rec_dose<-which.min(mapply(function(p,c) optimise(distance, c(0,target_pat), alpha=alpha, pred_x= p, pred_y=c, patient_target_DLT=target_pat, clinician_target_DLT=target_clin)$objective, p=p_p, c=p_c))
    #rec_dose<-min(which.min(abs(p_c-target_clin)), which.min(abs(p_p-target_pat)))
    ifelse(rec_dose>max(M[,2]), current_dose<-max(M[,2])+1, current_dose<-rec_dose)
    val<-t(sapply(rep(current_dose, times=no_enrolled), function (k) time_to_dlt(k,true_tox_c, true_tox_p, phi)))
    M<-truncated_matrix(no_enrolled, val, M, current_dose)
  }
   est_c<-optimise(like_c,c(0,10), u=u_skeleton, data_frame=M, maximum = TRUE)
   est_p<-optimise(like_p,c(0,10), v=v_skeleton, data_frame=M, maximum = TRUE)
   #dose-finding decision 
   p_c<-u_skeleton^est_c$maximum 
   p_p<-v_skeleton^est_p$maximum
   rec_dose<-which.min(mapply(function(p,c) optimise(distance, c(0,target_pat), alpha=alpha, pred_x= p, pred_y=c, patient_target_DLT=target_pat, clinician_target_DLT=target_clin)$objective, p=p_p, c=p_c))
  ifelse(rec_dose>max(M[,2]), current_dose<-max(M[,2])+1, current_dose<-rec_dose)
  unlink(file.path("C:/Users/ealger/AppData/Local/Temp", "Rtmp*"), recursive = T)
  optimal_dose<- min(which.min(abs(true_tox_c-target_clin)), which.min(abs(true_tox_p-target_pat)))
  mtd_assigned<-sum(M[,2]==optimal_dose)/sample
  overdosed<- sum(M[,2]>optimal_dose)/sample
  return(c(current_dose, mtd_assigned,overdosed))
}

