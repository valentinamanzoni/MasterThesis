

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

normboot.stpm2 <- function(fit,B){
  require(MASS)
  
  res <- mvrnorm(B, mu=fit@coef, Sigma=fit@vcov)
  res <- res[,order(colnames(res))] 
  return(res)
  
}

getknots <- function(fit){
  timeindex <- grep("bf",names(fit@lm$model))
  knots <- attributes(fit@lm$model[[timeindex]])$knots
  bknots <- attributes(fit@lm$model[[timeindex]])$Boundary.knots
  
  return(c(bknots[1],knots))
}


fix_coefs <- function(i,obj,tmat,basepar,transcovs){
  basepar[[i]] <- basepar[[i]][-length(basepar[[i]])]
  # print(basepar)
  # print(i)
  names_old <- colnames(obj[[i]])
  
  #identify to
  to_coef <- grep("to",substr(names_old,1,2))
  names_to <- substr(names_old[to_coef],3,10000)
  #identify from
  from_coef <- grep("from",substr(names_old,1,4))
  names_from <- substr(names_old[from_coef],5,10000)
  
  # we need to remove coefficients that are not of the right transition
  if(length(names_to)>0){
    trans <- colnames(tmat)[which(tmat == i, arr.ind = TRUE)[2]]
    index_to_keep <- to_coef[grepl(trans,names_to)]
    index_to_remove <- setdiff(to_coef,index_to_keep)
  }
  
  if(length(names_from)>0){
    trans <- colnames(tmat)[which(tmat == i, arr.ind = TRUE)[1]]
    index_to_keep <- from_coef[grepl(trans,names_from)]
    index_to_remove <- setdiff(from_coef,index_to_keep)
  }
  
  obj_new <- obj
  obj_new[[i]] <- obj[[i]][,-index_to_remove,drop=F]
  names_old <- colnames(obj_new[[i]])
  
  # identify baseline
  gamma0 <- grep("Intercept",names_old)
  others_gamma <- grep("bf",names_old)
  base_others_gamma <- others_gamma[1:max(basepar[[i]])]
  
  #rename baseline
  colnames(obj_new[[i]])[c(gamma0,base_others_gamma)] <- paste0("gamma",0:length(base_others_gamma))
  
  if(length(index_to_keep)>0){
    nocovs <- c(grep("gamma",colnames(obj_new[[i]])),grep("to",substr(colnames(obj_new[[i]]),1,2)),grep("from",substr(colnames(obj_new[[i]]),1,4)))
    covs <- colnames(obj_new[[i]])[-nocovs]
    #identify to
    to_coef <- grep("to",substr(colnames(obj_new[[i]]),1,2))
    names_to <- substr(colnames(obj_new[[i]])[to_coef],3,10000)
    #identify from
    from_coef <- grep("from",substr(colnames(obj_new[[i]]),1,4))
    names_from <- substr(colnames(obj_new[[i]])[from_coef],5,10000)
    
    #these are the coef relative to covariates that needs to be sum up
    names_to_covs <- lapply(covs,function(x)to_coef[grep(x,names_to)])
    names_from_covs <- lapply(covs,function(x)from_coef[grep(x,names_from)])
    
    if(length(to_coef)>0){
      # let's sum up the baseline 
      if(length(unlist(names_to_covs))>0){
        to_base <- to_coef[to_coef!=unlist(names_to_covs)]
      }else{
        
        to_base <- to_coef
      }
      
      
      base <- grep("gamma",colnames(obj_new[[i]]))
      if(length(to_base)!=length(base))errorCondition("Something is wrong in the summing up of the baseline!") 
      obj_new[[i]][,base]<- obj_new[[i]][,base]+obj_new[[i]][,to_base]
      
      
      for(x in 1:length(covs)){
        obj_new[[i]]<- fix_trans_covs(x,obj_new[[i]],covs,names_to_covs,transcovs)
      }
      obj_new[[i]] <-obj_new[[i]][,-c(to_base,unlist(names_to_covs))] 
    }
    
    if(length(from_coef)>0){
      # let's sum up the baseline 
      if(length(unlist(names_from_covs))>0){
        from_base <- from_coef[from_coef!=unlist(names_from_covs)]
      }else{
        
        from_base <- from_coef
      }
      base <- grep("gamma",colnames(obj_new[[i]]))
      if(length(from_base)!=length(base))errorCondition("Something is wrong in the summing up of the baseline!") 
      obj_new[[i]][,base]<- obj_new[[i]][,base]+obj_new[[i]][,from_base]
      for(x in 1:length(covs)){
        obj_new[[i]]<- fix_trans_covs(x,obj_new[[i]],covs,names_from_covs,transcovs)
      }
      obj_new[[i]] <-obj_new[[i]][,-c(from_base,unlist(names_from_covs)),drop=F] 
    }
  }
  if(is.null(dim(obj_new[[i]]))) {
    m <- matrix(obj_new[[i]],nrow=1)
    colnames(m) <- names(obj_new[[i]])
   
    return(m)
  }else
    return(obj_new[[i]])
}


fix_trans_covs <- function(x,f,covs,names,transcovs){
  index_covs <- unlist(names[[x]])
  if(is_empty(index_covs)){
    return(f)
  }else{
    if(!(x%in%transcovs)){
      f <- f[,-index_covs,drop=F]
    }else{
      index <- grep(covs[[x]],colnames(f))
      index <- index[index!=index_covs]
      f[,index]<- f[,index,drop=F]+f[,index_covs,drop=F]
    }
    
    return(f)
  }
  
  
}

bf <- function (x, knots){
  res <- basis(knots, x)[,-c(1,length(knots))]
  attributes(res)$knots <- knots[-c(1,length(knots))]
  attributes(res)$Boundary.knots <- knots[c(1,length(knots))]
  return(res)
}




prepare_models_list_boot <- function(model_mm,model_death,tmat,B){
  fits <- list(model_mm,
               model_death)
  fboot <- lapply(1:2,function(x) {
    n <- normboot.stpm2(fits[[x]],B=B)
    names <- colnames(n)
    n <- rbind(fits[[x]]@coef,n)
    colnames(n) <- names
    return(n)})
  
  
  death_trans <- na.omit(tmat[,ncol(tmat)])
  fit_list <- rep(list(fboot[[1]]), times = max(tmat,na.rm = T))
  fit_list[death_trans] <-rep(list(fboot[[2]]), times = length(death_trans))
  
  return(fit_list)
}


prepare_models_list_mle <- function(model_mm,model_death,tmat,B){
  fits <- list(model_mm,
               model_death)
  fboot <- lapply(1:2,function(x) {
    m <- matrix(fits[[x]]@coef,nrow=1)
    colnames(m) <- names(fits[[x]]@coef)
    return(m)
  })
  
  
  death_trans <- na.omit(tmat[,ncol(tmat)])
  fit_list <- rep(list(fboot[[1]]), times = max(tmat,na.rm = T))
  fit_list[death_trans] <-rep(list(fboot[[2]]), times = length(death_trans))
  
  return(fit_list)
}


prepare_models_list <- function(model_mm,model_death,tmat){
  death_trans <- na.omit(tmat[,ncol(tmat)])
  fit_list <- rep(list(model_mm), times = max(tmat,na.rm = T))
  fit_list[death_trans] <- model_death
  
  return(fit_list)
}


get_params <- function(model_mm,
                       model_death,tmat,B,transcovs){
  
  fboot <-prepare_models_list_boot(model_mm = model_mm,
                                   model_death = model_death,
                                   tmat=tmat,
                                   B=B)
  
  knots <- lapply(prepare_models_list(model_mm = model_mm,
                                      model_death = model_death,
                                      tmat=tmat),function(x)getknots(x))
  
  basepar <- lapply(knots, function(x)1:(length(x)))
  fits_fixed <- lapply(1:length(fboot),function(x)fix_coefs(x,fboot,tmat,basepar,transcovs))
  coefs_psa <-lapply(1:length(fboot) , function(x)coef_to_list_boot(fits_fixed [[x]] ,basepar[[x]]))
  transmod_params_psa <- define_params(coefs_psa,knots)
  return(list(transmod_params_psa=transmod_params_psa,fboot=fboot,fits_fixed=fits_fixed))
}


get_params_mle <- function(model_mm,
                           model_death,tmat,B,transcovs,n){
  
  fboot <-prepare_models_list_boot(model_mm = model_mm,
                                   model_death = model_death,
                                   tmat=tmat,
                                   B=B)
  
  knots <- lapply(prepare_models_list(model_mm = model_mm,
                                      model_death = model_death,
                                      tmat=tmat),function(x)getknots(x))
  
  basepar <- lapply(knots, function(x)1:(length(x)))
  fits_fixed <- lapply(1:length(fboot),function(x)fix_coefs(x,fboot,tmat,basepar,transcovs))
  coefs_psa <-lapply(1:length(fboot) , function(x)coef_to_list_boot(fits_fixed [[x]] ,basepar[[x]],n=n))
  transmod_params_psa <- define_params(coefs_psa,knots)
  return(list(transmod_params_psa=transmod_params_psa,fboot=fboot,fits_fixed=fits_fixed))
}







define_params <- function(coefs,knots) {
  
  coefage1<- matrix(1,
                    nrow=1,
                    ncol=1)
  colnames(coefage1) <- "Age_entry1"
  
  
  coefage2<- matrix(1,
                    nrow=1,
                    ncol=1)
  colnames(coefage2) <- "Age_entry2"
  
  
  coefage3<- matrix(1,
                    nrow=1,
                    ncol=1)
  colnames(coefage3) <- "Age_entry3"
  
  
  coefage4<- matrix(1,
                    nrow=1,
                    ncol=1)
  colnames(coefage4) <- "Age_entry4"
  
  coefage5<- matrix(1,
                    nrow=1,
                    ncol=1)
  colnames(coefage5) <- "Age_entry5"

  #  #initialize object
  transmod_params <- params_surv_list(
    # 0. from 0 to age entry
    params_surv(coefs = list(est = coefage1),
                dist = "fixed"),
    params_surv(coefs = list(est = coefage2),
                dist = "fixed"),
    params_surv(coefs = list(est = coefage3),
                dist = "fixed"),
    params_surv(coefs = list(est = coefage4),
                dist = "fixed"),
    params_surv(coefs = list(est = coefage5),
                dist = "fixed"))
  
  
  # add for all transitions
  for (trans in 1:length(coefs)){
    transmod_params[[trans+5]] <- params_surv(
      coefs = coefs[[trans]],
      dist = "survspline",
      aux = list(
        knots = knots[[trans]],
        scale = "log_cumhazard",
        timescale = "log"
      ))
  }
  
  return(transmod_params)
}

coef_to_list <- function(fit){
  
  
  
  coeflist <- lapply(fit$basepars, fix_par,fit=fit)
  names(coeflist) <- paste0("gamma",0:(length(coeflist)-1))
  return(coeflist)
}





fix_par <- function(index,fit){
  
  if(index==1){
    pars <- t(as.matrix(fit$coefficients[grepl("gamma0",names(fit$coefficients))|!grepl("gamma",names(fit$coefficients))  ]))
  }else{
    pars <- t(as.matrix(fit$coefficients[grepl(paste0("gamma",index-1),names(fit$coefficients))]))
  }
  
  if(ncol(pars)==1)colnames(pars) <- "intercept"
  if(ncol(pars)!=1){
    if(index==1){
      colnames(pars) <- c(
        "intercept",
        names(fit$coefficients)[!grepl("gamma",names(fit$coefficients))])
    }else{
      colnames(pars) <-c(
        "intercept",
        gsub('[()]',"",gsub(paste0("gamma",index-1),"",names(fit$coefficients)[grepl(paste0("gamma",index-1),names(fit$coefficients))][-1]))) 
    }
  }
  
  return(pars)
}



fix_par_boot <- function(index,coefboot){
  
  if(index==1){
    pars <- as.matrix(coefboot[,grepl("gamma0",colnames(coefboot))|!grepl("gamma",colnames(coefboot))  ])
  }else{
    pars <- as.matrix(coefboot[,grepl(paste0("gamma",index-1),colnames(coefboot))])
  }
  
  if(ncol(pars)==1)colnames(pars) <- "intercept"
  if(ncol(pars)!=1){
    if(index==1){
      colnames(pars) <- c(
        "intercept",
        colnames(coefboot)[!grepl("gamma",colnames(coefboot))])
    }else{
      col <- c("intercept",
               gsub('[()]',"",gsub(paste0("gamma",index-1),"",colnames(coefboot)[grepl(paste0("gamma",index-1),colnames(coefboot))][-1]))) 
      
      colnames(pars) <-col
    }
  }
  
  return(pars)
}




coef_to_list_boot <- function(fit,basepars,n){
  
  
  coeflist <- lapply(basepars, function(x){
    res <- fix_par_boot(x,coefboot=fit)
    res_sel <- res[1:n,,drop=F]
    colnames(res_sel) <- colnames(res)
    return(res_sel)
  }
  )
  names(coeflist) <- paste0("gamma",0:(length(coeflist)-1))
  return(coeflist)
}





sim_life <- function(patients,strategies,states,tmat,transmod_params_psa,B=1){
  
  gc()
  
  hesim_dat <- hesim_data(strategies = strategies,
                          patients = patients)
  
  transmod_data <-hesim::expand(hesim_dat, by = c("strategies", "patients"))
  
  tmat <- cbind(rep(NA,ncol(tmat)+1),rbind(c(1:(ncol(tmat)-1),NA),tmat+(ncol(tmat)-1)))
  colnames(tmat)[1] <- rownames(tmat)[1] <- "Pre"

  
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
    clock = "forward",
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
