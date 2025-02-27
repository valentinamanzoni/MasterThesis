

#data is in mstate format
get_states <- function(tmat){
  
  states <- rownames(tmat)
  
  dat_states <- data.table(state_id=1:length(states),
                           state_name=states)
  return(dat_states)
  
}


get_tmat <- function(data){
  tm <- t(as.matrix(table(data$from,data$to)))
  tm[tm==0] <- NA
  tm[!is.na(tm)] <- cumsum(!is.na(tm))[!is.na(tm)]
  tm <- t(tm)
  class(tm) <- "matrix"
  return(tm)
}






simulate_mstateB <- function(scenario_obj){
  
}


sim_life <- function(patients,strategies,states,tmat,transmod_params_psa,transition_types,B=1){
  
  gc()
  
  tmat <- cbind(rep(NA,ncol(tmat)+1),rbind(c(1:(ncol(tmat)-1),NA),tmat+(ncol(tmat)-1)))
  colnames(tmat)[1] <- rownames(tmat)[1] <- "Pre"
  

  hesim_dat <- hesim_data(strategies = strategies,
                          patients = patients)
  
  transmod_data <-hesim::expand(hesim_dat, by = c("strategies", "patients"))
  
  states <- get_states(tmat)
  # Utilities
  utility_tbl <- stateval_tbl(
    data.table(state_id=states$state_id[-length(states)],
               est=1),
    dist = "fixed"
  )
  
  # Model for transitions
  transmod <-create_IndivCtstmTrans(
    transmod_params_psa,
    input_data = transmod_data,
    trans_mat = tmat,
    clock = "forward" ,
    start_age = 0
  )

  # Model for utilities
  utilitymod <- create_StateVals(utility_tbl,
                                 n=B,
                                 hesim_data =hesim_dat)
  
  
  # MODEL
  ictsm <- IndivCtstm$new(trans_model = transmod,
                          utility_model = utilitymod)
  
  ictsm$sim_disease(
    max_age = 120,
    max_t=120,
    progress=F)

  return(ictsm$disprog_)
  
  
}

fix_subjects <- function(covs,f,nsim){
  
  dat <-
    as.data.frame(model.matrix(
      as.formula(paste0("~",f)),
      data = covs
    ))
  
  colnames(dat)[1] <-"intercept"
 
  
  return(dat)
}
