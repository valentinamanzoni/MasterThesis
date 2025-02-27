
sim_mm_traj <- function(n,N,scenario_obj,seed,scenario){
  
  require(tictoc)
  require(parallel)
  
  source(paste0("Simulation//fun//utils_microsim",scenario,".R"))
  

  set.seed(seed)
  nclust <- max(scenario_obj$pattern_obj$obj$predclass)
  tmat <- scenario_obj$tmat
  mstate_fits <- scenario_obj$trans_mod
  dis_fits <- scenario_obj$incid_dis
  dis_fits_rare <- scenario_obj$inc_dis_rare
 
  diseases <- get_diseases(scenario_obj)
  rare <- names(scenario_obj$inc_dis_rare)
  prev_diseases <- get_prevalence(scenario_obj)
  colnames(scenario_obj$prev_dis)[2:62] <- tolower(colnames(scenario_obj$prev_dis)[2:62])
 

###########################################################
######   1. SIMULATE LATENT MULTIMORBIDITY PATTERNS #######
###########################################################
  
  suppressMessages({
    suppressWarnings({
    # Patients
    covs <- data.frame(patient_id=1:(n*N),
                       id=rep(1:n,times=N),
                       educ_el=rbinom(n*N,1,prob=0.15),
                       dm_sex=rbinom(n*N,1,prob=0.45),
                       grp_id=rep(1:N,each=n),
                       Age_entry=pmin(rgamma(n*N, shape=0.9, rate = 0.15) + 60,96),
                       time_in_study= runif(n*N,0.5,20))#to modify according schema

    covs$start_age_cut <- cut(covs$Age_entry,breaks=c(0,66,72,78,81,110),include.lowest=T)
    print(prop.table(table(covs$start_age_cut)))
    print(summary(covs$Age_entry))
   
    covs$initial_state <- get_initial_mm_state(covs,scenario_obj$prev_initial_latent_state)
    
    covs$educ_el_mm=covs$educ_el
    covs$dm_sex_mm=covs$dm_sex
    
    f <- "patient_id+grp_id+educ_el+dm_sex+educ_el_mm+dm_sex_mm+Age_entry+initial_state+time_in_study+id"
    patients <- fix_subjects(covs,f = f,nsim =n*N)
    # Scenario
    strategies <- data.frame(strategy_id = 1)
  
    ###########################################################
    ######   2. SIMULATE LATENT MULTIMORBIDITY PATTERNS #######
    ###########################################################
    
    if (scenario=="B") {
      
      entry1 <- matrix(1,
                      nrow=1,
                      ncol=1)
      
      colnames(entry1) <- "Age_entry1"
      
      entry2 <- matrix(1,
                       nrow=1,
                       ncol=1)
      
      colnames(entry2) <- "Age_entry2"
      
      patients$Age_entry1 <- ifelse(patients$initial_state==1,patients$Age_entry,200)
      patients$Age_entry2 <- ifelse(patients$initial_state==2,patients$Age_entry,200)
      patients$Intercept <- 1
      colnames(patients)[ncol(patients)] <- "(Intercept)"
      colnames(patients)[colnames(patients)=="dm_sex"] <- "dm_sexFemales"
      
      params <- params_surv_list(params_surv(coefs = list(est=as.matrix(entry1)),dist = "fixed"),
                                 params_surv(coefs = list(est=as.matrix(entry2)),dist = "fixed"),
                                 create_params(mstate_fits[[1]] ,uncertainty = "none"),
                                 create_params(mstate_fits[[2]] ,uncertainty = "none"),
                                 create_params(mstate_fits[[3]] ,uncertainty = "none"))
      
  
      tic("MM pattern")

      mmtraj <- sim_life(patients = patients,
                         strategies = strategies,
                         states=c("Pre",states),
                         tmat=tmat,
                         transmod_params_psa=params)
      
      
    }else{
      params <- get_params_mle(model_mm = mstate_fits[[1]],
                               model_death = mstate_fits[[2]],
                               tmat = tmat,
                               B=2,
                               n=1,
                               transcovs=c(1,2,3,4,5,9,12,16,19))
      
      patients$Age_entry1 <- ifelse(patients$initial_state==1,patients$Age_entry,200)
      patients$Age_entry2 <- ifelse(patients$initial_state==2,patients$Age_entry,200)
      patients$Age_entry3 <- ifelse(patients$initial_state==3,patients$Age_entry,200)
      patients$Age_entry4 <- ifelse(patients$initial_state==4,patients$Age_entry,200)
      patients$Age_entry5 <- ifelse(patients$initial_state==5,patients$Age_entry,200)

      tic("MM pattern")
    
      mmtraj <- sim_life(patients = patients,
                         strategies = strategies,
                         states=states,
                         tmat=tmat,
                         transmod_params_psa=params$transmod_params_psa,
                         B=2)
      toc()
      
    }
    

    error <- mmtraj %>% filter(time_start>time_stop) %>% nrow()
    
    if(error>0) message("Problem with hesim simulation...")
    
    gc()
    
    mmtraj%<>% as.data.frame()%>% 
      left_join(covs) %>% 
      mutate(patient_id=id) %>% 
      dplyr::select(-id) %>% 
      filter(!(from==1)) %>% 
      rename(dataset_id=grp_id) %>% 
      mutate(from=from-1,
             to=to-1) %>% 
      rename(MP=from,
             Age=time_start) %>% 
      dplyr::group_by(dataset_id,patient_id) %>% 
      mutate(Age_death=max(time_stop[to==ncol(tmat)],na.rm = T),
             Age_death=na_if(Age_death,-Inf)) %>% 
      dplyr::select(dataset_id,patient_id,Age,MP,Age_death,Age_entry,educ_el,dm_sex,Age_entry,time_in_study)
      
    
    
    
    #########################################################
    ######   2. SIMULATE  PREVALENT DISEASES          #######
    #########################################################

    MM_baseline <- mmtraj %>% 
      arrange(dataset_id,patient_id) %>% 
      filter(Age==Age_entry) %>% 
      dplyr::select(patient_id,dataset_id,MP,Age) %>% 
      rename(base_MP=MP,
             Age_base_MP=Age)
    
    
    
    mmtraj%<>%
      left_join(MM_baseline) %>% 
      group_by(patient_id,dataset_id) %>% 
      mutate(MP=ifelse(Age>(time_in_study+Age_entry),lag(MP,default = 1),MP))%>% 
      mutate(Age=ifelse(Age>(time_in_study+Age_entry),time_in_study+Age_entry,Age)) %>% 
      distinct(patient_id,dataset_id,Age,.keep_all = T)
    
    exit <- mmtraj %>%
      #filter((Age_entry+time_in_study)>Age_death) %>%
      group_by(dataset_id,patient_id) %>%
      filter(row_number()==n()) %>%
      mutate( Age=pmin(Age_death,Age_entry+time_in_study,na.rm = T))

    mmtraj %<>%
      bind_rows(mmtraj,exit) %>%
      arrange(dataset_id,patient_id,Age) %>% 
      distinct(dataset_id,patient_id,Age,.keep_all=T)
   
    
    if(length(diseases)!=length(prev_diseases)) 
      errorCondition("Check the list of chronic conditions! It does not match with the prevalences...")
    tic("Simulated prevalences")
    
    
    
    # Convert base_MM to data.table for faster operations
    MM_baseline <- as.data.table(MM_baseline)
    MM_baseline %<>% group_by(patient_id,dataset_id) %>% mutate(id=cur_group_id()) 
    MM_baseline %<>% mutate(Age_cut=cut(Age_base_MP,breaks = c(0,65,70,75,80,85,90,95,110),
                                        include.lowest = T)) 
    
    MM_baseline %<>% arrange(dataset_id,patient_id)
    
    
    simulated_prev_diseases <- matrix(simulate_prevalent_diseases_vectorized(diseases,
                                                    prev_diseases,MM_baseline),nrow = n*N,ncol=length(diseases),byrow = F)
    
    colnames(simulated_prev_diseases) <- tolower(diseases)
    
    
    simulated_prev_diseases <- cbind(simulated_prev_diseases,
                                     data.frame(dataset_id=MM_baseline$dataset_id,
                                                patient_id=MM_baseline$patient_id))
    
    
    simulated_prev_diseases_rare <- matrix(simulate_prevalent_diseases_rare_vectorized(rare,
                                                                                       scenario_obj$prev_dis,
                                                                                       MM_baseline),nrow = n*N,ncol=length(rare),byrow = F)
      
      
    
    colnames(simulated_prev_diseases_rare) <- rare
    
    simulated_prev_diseases_rare <- cbind(simulated_prev_diseases_rare,
                                     data.frame(dataset_id=MM_baseline$dataset_id,
                                                patient_id=MM_baseline$patient_id))
    
    simulated_prev_diseases %<>% full_join(simulated_prev_diseases_rare)
    simulated_prev_diseases %<>% dplyr::select(dataset_id,patient_id,everything()) 
    toc()
    #########################################
    ######   3. REMOVE NO MM          #######
    #########################################
    
    # simulare anche le malattie rare e poi eliminare i soggetti non MM
    
    ndis <- apply(simulated_prev_diseases[,3:61],1,sum)
    dat_n_dis <- data.frame(dataset_id=MM_baseline$dataset_id,
                            patient_id=MM_baseline$patient_id,
                            ndis=ndis)
    simulated_prev_diseases %<>% left_join(dat_n_dis) %>% 
      arrange(dataset_id,patient_id) %>% filter(ndis>1)
   
    
    mmtraj %<>%left_join(dat_n_dis) %>% 
      left_join(simulated_prev_diseases) %>% 
      filter(ndis>1) %>% 
      arrange(dataset_id,patient_id,Age) %>% 
      group_by(patient_id,dataset_id) %>% 
      mutate(lag_age=lag(Age))
    
    mmtraj %>% distinct(patient_id, dataset_id) %>% nrow() %>% print()
    for (i in 1:N){
      print(paste(n-length(unique(mmtraj$patient_id[mmtraj$dataset_id==i])),"non-mm subjects removed in dataset",i,"..."))
    }
    
   print(nrow(simulated_prev_diseases))
    
    #######################################################
    ######   4. SIMULATE INCIDENT DISEASES          #######
    #######################################################
    
    tic("Simulated incident diseases...")
    
    simulated_inc_diseases <- matrix(simulate_incident_diseases_vectorized(diseases,
                                                                           mmtraj=mmtraj,
                                                                           dis_fits=dis_fits,
                                                                           prev_diseases = prev_diseases),
                                     ncol=length(diseases),
                                     byrow=F,nrow=nrow(simulated_prev_diseases))
    
    colnames(simulated_inc_diseases) <- paste0("Age_",tolower(diseases))
    dat_simulated_inc_diseases <-  cbind(as.data.frame(simulated_inc_diseases),
                                        data.frame(dataset_id=simulated_prev_diseases$dataset_id,
                                                   patient_id=simulated_prev_diseases$patient_id))
    
    simulated_inc_diseases_rare <- matrix(simulate_incident_diseases_rare_vectorized(rare[rare!="chromosomal_abnormalities"],mmtraj=mmtraj,dis_fits=dis_fits_rare),ncol=length(rare),byrow=F,nrow=nrow(simulated_prev_diseases))
    
    print(paste0("Age_",tolower(rare)[tolower(rare)!="chromosomal_abnormalities"]))
    colnames(simulated_inc_diseases_rare) <- paste0("Age_",tolower(rare)[rare!="chromosomal_abnormalities"])
    
    dat_simulated_inc_diseases_rare <- cbind(simulated_inc_diseases_rare,
                                             data.frame(dataset_id=simulated_prev_diseases$dataset_id,
                                                        patient_id=simulated_prev_diseases$patient_id))
    dat_simulated_inc_diseases %<>% full_join(dat_simulated_inc_diseases_rare)
    gc()
    toc()
    ###############################################################################
    ######   5. COMBINE DATA  FOR EXACTLY OBSERVED DATA                     #######
    ###############################################################################
    
    
    dat_res <- mmtraj %>% 
      distinct(dataset_id, patient_id, Age_entry,Age_death,.keep_all = T) %>% 
      dplyr::select(-Age,-MP,-lag_age) %>% 
      arrange(dataset_id, patient_id) 
    
    
    col_dis_name <- tolower(diseases)
    col_dis_name <- c(col_dis_name,rare[rare!="chromosomal_abnormalities"])
    #print(col_dis_name)
    dat_res[,col_dis_name] <- dat_res[,col_dis_name]*dat_res$Age_entry
    colnames(dat_res)[colnames(dat_res) %in% col_dis_name] <- paste0("Age_", col_dis_name)
    col_dis_name <- paste0("Age_", col_dis_name)
    dat_res %<>%
      mutate_at(col_dis_name,na_if,0) %>%
      left_join(dat_simulated_inc_diseases,by = c("patient_id","dataset_id"), suffix = c("_P", "_I"))%>%
      mutate(across(ends_with("_P"), ~ coalesce(.x, get(sub("_P", "_I", cur_column()))))) %>%
      dplyr::select( -ends_with("_I")) %>%
      rename_with(~ sub("_P", "", .x))

    
    #############
    # dat_res[,9:67] <- dat_res[,9:67]*dat_res$Age_entry
    # colnames(dat_res)[9:67] <- paste0("Age_",colnames(dat_res[,9:67]))
    # 
    # dat_res %<>%
    #   mutate_at(9:67,na_if,0) %>%
    #   left_join(dat_simulated_inc_diseases,by = c("patient_id","dataset_id"), suffix = c("_P", "_I"))%>%
    #   mutate(across(ends_with("_P"), ~ coalesce(.x, get(sub("_P", "_I", cur_column()))))) %>%
    #   dplyr::select( -ends_with("_I")) %>%
    #   rename_with(~ sub("_P", "", .x))

    }) })
  
  # simulare con unif centrata su vent'anni un tempo nello studio e poi rimuovere tutte le malattie dopo
  # togliere morte (lasciare ma ggiungere binaria 0-1 per morte) per persone che stanno nello studio meno di quando viene osservata morte
  dat_res$Age_exit <- pmin(dat_res$Age_death, dat_res$time_in_study+dat_res$Age_entry )
  MM_baseline %<>% dplyr::select(patient_id,dataset_id,base_MP)%<>% rename(MP_base=base_MP)
  dat_res %<>%left_join(MM_baseline, by = c("patient_id", "dataset_id"))
  #MM_baseline_selected <- MM_baseline %>% dplyr::select(patient_id, dataset_id, base_MP)s
  #MM_baseline_renamed <- MM_baseline_selected %>% rename(MP_base = base_MP)
  #dat_res <- dat_res %>% left_join(MM_baseline_renamed, by = c("patient_id", "dataset_id"))
  gc()
  return(list(dat_res, mmtraj))
}


simulate_prevalent_diseases <- function(disease,prev_diseases,base_MM){
  print("Real prevalence")
  print(prev_diseases[[disease]])
  probs <- prev_diseases[[disease]][base_MM$base_MP]
  dis <-  rbinom(nrow(base_MM), 1, prob = probs)
  print("Simulated")
  print(tapply(dis,base_MM$base_MP,function(x)prop.table(table(x))[2]))
  print(disease)
  return(dis)
}

simulate_prevalent_diseases_vectorized <- Vectorize(simulate_prevalent_diseases,"disease")



simulate_prevalent_diseases_rare <- function(disease,prev_diseases,base_MM){
  
  probs <- prev_diseases[[disease]][base_MM$Age_cut]
  dis <- rbinom(nrow(base_MM),1,prob=probs)
  print(disease)
  return(dis)
}

simulate_prevalent_diseases_rare_vectorized <- Vectorize(simulate_prevalent_diseases_rare,"disease")


rtrunc_v <- Vectorize(rtrunc,c("a","b"))
simulate_incident_diseases <- function(disease,mmtraj,dis_fits,prev_diseases){
  
  
  
  mmtraj %<>% drop_na(lag_age) 
  
  require(flexsurv)
  require(truncdist)
  
  idnoprev <- which(mmtraj[tolower(disease)]==0)
  
  print("Real prevalence")
  print(prev_diseases[[disease]])
  
  probs <- prev_diseases[[disease]][mmtraj$MP]
  mmtraj$dis <-  rbinom(nrow(mmtraj), 1, prob = probs)
  
  print("Simulated")
  print(tapply(mmtraj$dis,mmtraj$MP,function(x)prop.table(table(x))[2]))
  
  lower <- ifelse(mmtraj$Age == mmtraj$Age_death | mmtraj$Age ==(mmtraj$Age_entry + mmtraj$time_in_study),(mmtraj$lag_age-60)/62 ,pmax((mmtraj$lag_age-60)/62,(mmtraj$Age-63)/62))
  
  upper <- (mmtraj$Age-60)/62
  
  mmtraj$age_dis <- mapply(rtrunc,a=lower,b=upper,MoreArgs=list(n=1,
                                                spec="beta",
                               shape1=dis_fits[1,tolower(disease)],
                               shape2=dis_fits[2,tolower(disease)]))*62+60
  dat_dis <- mmtraj %>% 
    ungroup() %>%
    # group_by(dataset_id,patient_id) %>% 
    mutate(age_dis=ifelse(dis==1,age_dis,NA_integer_)) %>%
    #summarize(age_dis=min(age_dis,na.rm = T)) %>%
    #mutate(age_dis=ifelse(age_dis==Inf,NA_integer_,age_dis)) %>%
    dplyr::mutate(age_dis=ifelse(row_number()%in%idnoprev,age_dis,NA_integer_)) %>% 
    arrange(dataset_id,patient_id, age_dis) %>% 
    distinct(dataset_id,patient_id, .keep_all = TRUE)
  
  # dat_dis %>% filter(dataset_id==1 & patient_id==32) %>% dplyr::select(Age, age_dis) %>% print()
  # if (disease=="Dementia"){
  #   View(dat_dis)
  # }
  print(nrow(dat_dis))

    
  dat_dis %>% filter(dataset_id==1 & patient_id==32) %>% dplyr::select(age_dis) %>% print()
  dat_dis %<>%pull(age_dis)

  print(disease)
  return(dat_dis)
  
}


simulate_incident_diseases_vectorized <- Vectorize(simulate_incident_diseases,"disease")
simulate_incident_diseases_rare <- function(disease,mmtraj,dis_fits){
  require(flexsurv)
  
  
  dat_dis <- mmtraj %>% 
    filter(.data[[disease]]==0) %>% 
    group_by(patient_id,dataset_id) %>% 
    summarize(tstart=min(Age),
              tstop=min(Age_death),
              id=cur_group_id())
  
  sim_dis <- simulate(dis_fits[[casefold(disease)]],
                      nsim=1,
                      newdat=dat_dis,
                      start=dat_dis$tstart,
                      censtime=dat_dis$tstop,
                      tydy=T)
  
  dat_dis$dis_age <- ifelse(sim_dis[,2]==1,sim_dis[,1],NA_real_)
  
  dat_dis %<>% dplyr::select(patient_id,dataset_id,dis_age) 
  colnames(dat_dis)[ncol(dat_dis)] <- paste0("Age_",disease)
  
  
  dat_dis_full <- mmtraj %>% 
    distinct(patient_id,dataset_id)
  
  dat_dis %<>% full_join(dat_dis_full) %>% arrange(dataset_id,patient_id) %>% pull(contains(paste0("Age_",disease)))
  
  print(disease)
  
  return(dat_dis)
  
}


simulate_incident_diseases_rare_vectorized <- Vectorize(simulate_incident_diseases_rare,"disease")



get_diseases <- function(scenario_obj){
  names(scenario_obj$pattern_obj$obj$probs) 
}

get_prevalence <- function(scenario_obj){
  lapply(scenario_obj$pattern_obj$obj$probs,function(x){
    y <- as.numeric(x[,2])
  
    return(y[scenario_obj$right_order])})
}


get_initial_mm_state <- function(covs,prev_initial_latent_state){
  state <- do.call(c,lapply(1:nrow(covs),function(x)which(rmultinom(n=1,size=1,prob=prev_initial_latent_state$p[prev_initial_latent_state$Age_cut==covs$start_age_cut[x]])==1)))
 print(length(state))
 return(state)
  }




