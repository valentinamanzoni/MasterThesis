library(msm)
library(tictoc)


apply_msm <-function(pop_ms,misc,sim_obj, result_folder){
if (scenario %in% c("A", "Av2")){
  # 5 state + death
  q_matrix_1 <- rbind(c(-0.5,0.1,0.1,0.1,0.1,0.1), c(0.1,-0.6,0.1,0.1,0.1,0.1),
                    c(0,0.1,-0.4,0.1,0.1,0.1), c(0,0.1,0.1,-0.4,0.1,0.1),
                    c(0,0.1,0.1,0.1,-0.4,0.1), c(0,0,0,0,0,0))

  q_matrix <- sim_obj$tmat
  q_matrix[is.na(q_matrix)] <-0
  q0 <- crudeinits.msm(MP ~ age, subject=patient_id, data=pop_ms, qmatrix= q_matrix)
  q0_1 <- crudeinits.msm(MP ~ age, subject=patient_id, data=pop_ms, qmatrix= q_matrix_1)
  tic("msm age")
  model_age <- msm(MP ~ age, subject=patient_id, data=pop_ms, qmatrix= q0_1, deathexact = 6,
                  covariates = ~ age,  control = list(fnscale = 10000))
  t<-toc()
  model_obj_a <- list(
    model = model_age,
    time = t$toc - t$tic
  )
  
  save(model_obj_a, file =paste0(result_folder,"/model_age_scenario_",scenario,"_", nsim,"_",unique(pop_ms$dataset_id),".RData"))
  
  qmatrix_init <- model_age$Qmatrices$baseline 
  covinits_init <-  as.list(model_age$Qmatrices$age[model_age$Qmatrices$age > 0])

  # model with misclassification and no age covariate
  n <- dim(misc)[1] -1
  # number of non empty element in q - diagonal
  n1<- sum(q0 > 0)
  n2<- n1+ (n*n)-n
  tic("msm misc")
  model_misc <-msm(MP ~ age, subject=patient_id, data=pop_ms, qmatrix= q0, 
                   deathexact = 6,
                   ematrix=as.matrix(misc), 
                   initprobs = c(0.6,0.1,0.1,0.1,0.1,0),
                   fixedpars= c((n1+1):n2), # da modificare 
                   control = list(fnscale = 20000))
  t<- toc()
  model_obj_m <- list(
    model = model_misc,
    time = t$toc - t$tic
  )
  save(model_obj_m, file= paste0(result_folder,"/model_misc_scenario_",scenario,"_", nsim,"_",unique(pop_ms$dataset_id),".RData"))
  # model with misclassification and age covariate
  n1<- sum(q0 > 0) + sum(covinits_init!=0)
  n2<- n1+ (n*n)-n
  tic("msm age + misc")
  model_age_misc <- msm(MP ~ age, subject=patient_id, data=pop_ms, qmatrix= qmatrix_init, 
      covinits = covinits_init, # passare lista valori age a convergenza
      deathexact = 6,
      ematrix=misc, 
      initprobs = c(0.6,0.1,0.1,0.1,0.1,0),
      covariates = ~ age,
      fixedpars= c((n1+1):n2), # da modificare 
      control = list(fnscale = 20000))
  t<- toc()
  model_obj_m_a <- list(
    model = model_age_misc,
    time = t$toc - t$tic
  )
  save(model_obj_m_a, file =paste0(result_folder,"/model_age_misc_scenario_",scenario,"_", nsim,"_",unique(pop_ms$dataset_id),".RData"))
  
} 
  else if (scenario=="B"){
  # 2 state + death 
  q_matrix <- rbind(c(0,1,1),
                    c(0,0,1),
                    c(0,0,0))
  ### crude inizialization
  q0 <- crudeinits.msm(MP ~ age, subject=patient_id, data=pop_ms, qmatrix= q_matrix)
  q0
  
  tic("msm age + cov")
  model_age_cov <- msm(MP ~ age, subject=patient_id, data=pop_ms, qmatrix= q0, deathexact = 3,
                       covariates = ~ age + educ_el + dm_sex,  control = list(fnscale = 10000))
  # covinits_init <- as.list(model_age_cov$Qmatrices$age[model_age_cov$Qmatrices$age > 0])
  t<-toc()
  model_obj_a <- list(
    model = model_age_cov,
    time = t$toc - t$tic
  )
  
  save(model_obj_a, file =paste0(result_folder,"/model_age_scenario_",scenario,"_", nsim,"_",unique(pop_ms$dataset_id),".RData"))
  
  
  ### misc + age cov + 2 covariate
  n <- dim(misc)[1] -1
  # number of non empty element in q - diagonal
  n1<- sum(q0 > 0)
  #n2<- n1+ (n*n)-n
  n2 <- length(model_age_cov$estimates)
  qmatrix_init <- model_age_cov$Qmatrices$baseline 
  # model_age_cov$Qmatrices$baseline 
  covinits_init <-model_age_cov$estimates
  covinits_init <-as.list(covinits_init[(n1+1):n2])
  n1<- sum(q0 > 0) + sum(covinits_init!=0)
  n2<- n1+ (n*n)-n
  tic("msm age + cov misc")
  model_age_misc <- msm(MP ~ age, subject=patient_id, data=pop_ms, qmatrix= qmatrix_init, 
                        covinits = covinits_init, # passare lista valori age a convergenza
                        deathexact = 3,
                        initprobs = c(0.8,0.2,0),
                        covariates = ~ age + educ_el + dm_sex,
                        ematrix=misc, 
                        fixedpars= c((n1+1):n2), # da modificare 
                        control = list(fnscale = 20000))
  t<- toc()
  model_obj_m_a <- list(
    model = model_age_misc,
    time = t$toc - t$tic
  )
  save(model_obj_m_a, file =paste0(result_folder,"/model_age_misc_scenario_",scenario,"_", nsim,"_",unique(pop_ms$dataset_id),".RData"))
  
  }
  else
    stop("Insert a valid design type")

  print("msm executed")
}