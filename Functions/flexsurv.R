library(flexsurv)
library(mstate)
library(msm) 




apply_flexsurv <- function(pop_ms, scenario_obj, result_folder){
  
  q_matrix<-scenario_obj$tmat
  q_matrix[is.na(q_matrix)] <-0
  data <- msm2Surv(data = pop_ms, subject = "patient_id", time= "age", state="MP", Q = q_matrix)
  # transizioni non permesse da q_matrix finiscono in 0
  #data[which(data$trans==0),]
  
  if (scenario =="B"){
    n_trans <-length(unique(data$trans))
  }else {
    # to account for "impossible" transitions
    n_trans <- length(unique(data$trans))-1
  }
  trans_model <- list()
  tic("flexsurv")
  for (i in 1:(n_trans)){
    model_0 <- flexsurvreg(formula = Surv(Tstart,Tstop, status) ~ educ_el + dm_sex, subset=(trans==i),
                           data = data, dist = "gompertz")
    trans_model[[i]]<- model_0
    names(trans_model)[i] <- paste(unique(data$from[data$trans==i]), unique(data$to[data$trans==i]), sep ="->")
  }
  t<-toc()
  # create object with the transition correspondence, the trans_model and the q_matrix
  models_object <- list(
    transmodel = trans_model,
    qmatrix = q_matrix,
    time = t$toc - t$tic
  )
  save(models_object, file = paste0(result_folder,"/flexsurv_scenario_",scenario,"_", nsim,"_", 
                                                             unique(pop_ms$dataset_id),".RData"))
  print("flexsurv executed")
  }

apply_flexsurv_base <- function(pop_ms, scenario_obj, result_folder){
  
  q_matrix<-scenario_obj$tmat
  q_matrix[is.na(q_matrix)] <-0
  data <- msm2Surv(data = pop_ms, subject = "patient_id", time= "age", state="MP_sim", Q = q_matrix)
  # transizioni non permesse da q_matrix finiscono in 0
  #data[which(data$trans==0),]
  
  if (scenario =="B"){
    n_trans <-length(unique(data$trans))
  }else {
    # to account for "impossible" transitions
    n_trans <- length(unique(data$trans))-1
  }
  trans_model <- list()
  tic("flexsurv")
  for (i in 1:(n_trans)){
    model_0 <- flexsurvreg(formula = Surv(Tstart,Tstop, status) ~ educ_el + dm_sex, subset=(trans==i),
                           data = data, dist = "gompertz")
    trans_model[[i]]<- model_0
    names(trans_model)[i] <- paste(unique(data$from[data$trans==i]), unique(data$to[data$trans==i]), sep ="->")
  }
  t<-toc()
  # create object with the transition correspondence, the trans_model and the q_matrix
  models_object <- list(
    transmodel = trans_model,
    qmatrix = q_matrix,
    time = t$toc - t$tic
  )
  save(models_object, file = paste0(result_folder,"/flexsurv_base_scenario_",scenario,"_", nsim,"_", 
                                    unique(pop_ms$dataset_id),".RData"))
  print("flexsurv executed")
}

apply_flexsurv_imp <- function(pop_ms, scenario_obj, result_folder){
  
  q_matrix<-scenario_obj$tmat
  q_matrix[is.na(q_matrix)] <-0
  data <- msm2Surv(data = pop_ms, subject = "patient_id", time= "age", state="MP", Q = q_matrix)
  
  if (scenario =="B"){
    n_trans <-length(unique(data$trans))
  }else {
    # to account for "impossible" transitions
    n_trans <- length(unique(data$trans))-1
  }
  trans_model <- list()
  tic("flexsurv")
  for (i in 1:(n_trans)){
    model_0 <- flexsurvreg(formula = Surv(Tstart,Tstop, status) ~ educ_el + dm_sex, subset=(trans==i),
                           data = data, dist = "gompertz")
    trans_model[[i]]<- model_0
    names(trans_model)[i] <- paste(unique(data$from[data$trans==i]), unique(data$to[data$trans==i]), sep ="->")
  }
  t<-toc()
  # create object with the transition correspondence, the trans_model and the q_matrix
  models_object <- list(
    transmodel = trans_model,
    qmatrix = q_matrix,
    time = t$toc - t$tic
  )
  save(models_object, file = paste0(result_folder,"/flexsurv_imp_scenario_",scenario,"_", nsim,"_", 
                                    unique(pop_ms$dataset_id),"_",unique(pop_ms$imp_nr),".RData"))
  print("flexsurv executed")
}

