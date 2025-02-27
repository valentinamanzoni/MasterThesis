library(msm)
library(flexsurv)
library(nhm)

##
### calculate rate and shape from msm:

hazards_mine <- function(x, b.covariates, no.years, trans = NULL,SE=FALSE, CI = FALSE,
                         age.shift = 0) {

  model <- x  # Fitted MSM model

  # Handle transition matrix (convert vector to matrix if needed)
  if (is.vector(trans)) {
    trans <- matrix(trans, 1, 2)
  }

  # Create the transition matrix if not provided
  if (!is.null(trans)) {
    ntrans <- nrow(trans)  # Number of transitions
  } else {
    Q <- model$qmodel$imatrix  # Transition rate matrix
    ntrans <- sum(Q)  # Count the number of transitions
    trans <- matrix(NA, ntrans, 2)  # Initialize matrix for transitions
    index <- 1
    for (i in 1:nrow(Q)) {
      for (j in 1:nrow(Q)) {
        if (Q[i, j]) {
          trans[index, ] <- c(i, j)  # Add transitions to matrix
          index <- index + 1
        }
      }
    }
  }

  # Number of states in the model
  nstates <- nrow(model$Qmatrices$baseline)

  # Create a null transition matrix with zeros on the diagonal
  Q.null <- matrix(as.numeric(model$Qmatrices$baseline != 0), nstates, nstates)
  diag(Q.null) <- 0
  ntrans <- sum(Q.null)  # Count the number of transitions

  # Number of covariates (including age and other covariates)
  ncovs <- 1 + length(b.covariates)
  nbeta <- max(which(names(model$estimates) == "qcov"))

  if (nbeta != (ntrans * ncovs)) {
    stop("\nAll covariates in msm model should be specified.\n\n")
  }

  # Extract the age variable from the covariates
  b.age <- b.covariates$age

  # Create a sequence of ages for computation
  age.grid <- seq(b.age, b.age + no.years, by = 1/12)

  # Initialize an empty list to store the hazards
  hazards <- list()

  # Loop over each transition to calculate the hazards
  for (i in 1:ntrans) {
    haz <- rep(NA, length(age.grid))
    hazLB <- rep(NA, length(age.grid))
    hazUB <- rep(NA, length(age.grid))
    hazSE <- rep(NA, length(age.grid))

    # For each age, calculate the hazard for the current transition
    for (j in 1:length(age.grid)) {
      if (ncovs == 2) {
        covariates <- list(age = age.grid[j])
      } else {
        covariates <- c(age = age.grid[j], b.covariates[2:length(b.covariates)])
      }

      # Compute the hazard for the current transition (without confidence intervals)
      haz[j] <- qmatrix.msm(model, covariates = covariates, ci = "none")[trans[i, 1], trans[i, 2]]

      if (CI || SE) {
        # Compute confidence intervals if requested
        Q.CI <- qmatrix.msm(model, covariates = covariates, ci = "delta", cores=4)
        hazLB[j] <- Q.CI$L[trans[i, 1], trans[i, 2]]
        hazUB[j] <- Q.CI$U[trans[i, 1], trans[i, 2]]
        hazSE[j] <- Q.CI$SE[trans[i, 1], trans[i, 2]]
      }
    }

    # Store the hazards (and optionally the confidence intervals)
    hazard_data <- list(hazard = haz)
    if (CI) {
      hazard_data$hazardLB <- hazLB
      hazard_data$hazardUB <- hazUB
    }
    if (SE){
      hazard_data$hazardSE <-hazSE
    }

    # Name the transition for easier identification later
    hazard_name <- paste("Transition", trans[i, 1], "to", trans[i, 2])
    hazards[[hazard_name]] <- hazard_data
  }

  # Return the computed hazards (and confidence intervals if requested)
  return(hazards)
}
# 
# params_msm_age <- matrix(model.msm_age$estimates[4:12], nrow = 3, ncol = 3) #1:3 rate 12:15 age
# params_msm_age <- cbind(params_msm_age, exp(params_msm_age[,1]), exp(params_msm_age[,2]), exp(params_msm_age[,3]))
# colnames(params_msm_age) <- colnames(ground_truth_params)[3:8]
# 
# min_age <- min(model.msm_age$data[[1]]$age)
# max_age <- max(model.msm_age$data[[1]]$age)
# haz <- hazards_mine(model.msm_age, b.covariates = list(age = 0, cov1 = 0, cov2 = 0, cov3 = 0), no.years = 40)
# 
# # Assuming this hazards come from fitting a gompertz model I wanna retrieve for each transition
# # shape and rate value, parameters of the distribution
# # h(t)=rate*exp(shape*t)
# # log(h(t))= log(rate) + shape*t
# # can be seen as y(t)= a+b*t
# shape <- numeric()
# rate <- numeric()
# 
# for (i in 1:3){
#   age_grid <- seq(min_age, max_age , length.out = length(unlist(haz[[1]])))
#   y <- log(as.numeric(unlist(haz[[i]])))
#   reg_model <- lm(y ~ age_grid)
#   rate[[i]] <- reg_model$coefficients[1]
#   shape[[i]]<- reg_model$coefficients[2]
# }
# 
# params_msm_age <- cbind(rate,shape,params_msm_age)
# bias_msm_age <- compute_bias(params_msm_age, ground_truth_params)


get_estimates <- function(model, model_name, i,j=0){
  # 
  rate<- numeric()
  rate_se <-numeric()
  shape <- numeric()
  shape_se <- numeric()
  trans <-numeric()
  beta_educ <- numeric()
  beta_sex <-numeric()
  beta_educ_se <- numeric()
  beta_sex_se <- numeric()

  if(model_name %in% c("flexsurv", "flexsurv_base", "flexsurv_imp"))
  {   
    endings <- names(model$transmodel)
    for (ending in endings) {
      current_model <- model$transmodel[[ending]]
      trans<- c(trans, ending)
      rate <-c(rate, current_model$res.t["rate", "est"])
      rate_se <- c(rate_se, current_model$res.t["rate", "se"])
      shape <- c(shape, current_model$res.t["shape", "est"])
      shape_se <- c(shape_se, current_model$res.t["shape", "se"])
      beta_educ <- c(beta_educ, current_model$res.t["educ_el", "est"])
      beta_educ_se <- c(beta_educ_se, current_model$res.t["educ_el", "se"])
      beta_sex <- c(beta_sex, current_model$res.t["dm_sex", "est"])
      beta_sex_se <- c(beta_sex_se, current_model$res.t["dm_sex", "se"])
    }
    if(model_name=="flexsurv_base"){
      model_name <- "flexsurv_true"
    }
  } else if (model_name %in% c("model_age", "model_age_misc")){
    # simil gompertz
    
    haz <- hazards_mine(model$model,SE=TRUE, b.covariates = list(age = 60, educ_el = 0, dm_sex = 0), no.years = 40)
    #haz <- hazards_mine(model$model, b.covariates = list(age = 0, educ_el = 0, dm_sex = 0), no.years = 40)
    
    # Assuming this hazards come from fitting a gompertz model I wanna retrieve for each transition
    # shape and rate value, parameters of the distribution
    # h(t)=rate*exp(shape*t)
    # log(h(t))= log(rate) + shape*t
    # can be seen as y(t)= a+b*t
    
    # >0 element in Qmatrices
    ntrans <-  sum(model$model$qmodel$imatrix)
    trans <- rename_trans(names(haz))
    min_age <- min(model$model$data[[1]]$age)
    max_age <- max(model$model$data[[1]]$age)
    
    # SE for rate and shape are computed with bootstrapping
    n_boot <- 1000
    for (i in 1:ntrans){
      age_grid <- seq(min_age, max_age , length.out = length(unlist(haz[[1]]$hazard)))
      y <- log(as.numeric(unlist(haz[[i]]$hazard)))
      reg_model <- lm(y ~ age_grid)
      rate <- c(rate,reg_model$coefficients[1])
      shape<- c(shape, reg_model$coefficients[2])
      haz_se <- unlist(haz[[i]]$hazardSE)
      
      bootstrap_results <- replicate(n_boot, {
        y_noisy <- log(as.numeric(unlist(haz[[i]]$hazard))) + rnorm(length(y), 0, haz_se)
        indices <- sample(1:length(y), replace = TRUE)
        y_boot <- y_noisy[indices]
        age_boot <- age_grid[indices]
      
        boot_model <- lm(y_boot ~ age_boot)
        boot_model$coefficients
      })
      
      rate_se <- c(rate_se, sd(bootstrap_results[1, ]))
      shape_se <- c(shape_se, sd(bootstrap_results[2, ]))

      #rate_se <- c(rate_se, )
      #shape_se <- rep(0,length(trans))
    }
    std_errors <- sqrt(diag(model$model$covmat))
    beta_educ <- model$model$estimates[((ntrans*2)+1):(3*ntrans)]
    beta_educ_se <- std_errors[((ntrans*2)+1):(3*ntrans)]
    beta_sex <- model$model$estimates[((ntrans*3)+1):(4*ntrans)]
    beta_sex_se <- std_errors[((ntrans*3)+1):(4*ntrans)]
    if(model_name=="model_age"){
      model_name <- "msm_base"
    } else if (model_name=="model_age_misc"){
      model_name <- "msm_misc"
    }

  } else if (model_name %in% c("nhm_misc", "nhm_base")){
    np <- model$model$nstate
    trans <- model$model$parnames[1:np]
    cov_matrix <- solve(model$model$hess)
    std_errors <- sqrt(diag(cov_matrix))
    trans <- gsub("^Base: ", "", trans)
    rate <- model$model$par[1:np]
    rate_se <- std_errors[1:np]
    shape <- model$model$par[(np+1):(2*np)]
    shape_se <- std_errors[(np+1):(2*np)]
    beta_educ <- model$model$par[((np*2)+1):(3*np)]
    beta_educ_se <- std_errors[((np*2)+1):(3*np)]
    beta_sex <- model$model$par[((np*3)+1):(4*np)]
    beta_sex_se <- std_errors[((np*3)+1):(4*np)]
  }
  if (model_name=="flexsurv_imp"){
    est_obj <- data.frame(
      model = rep(model_name, length(trans)),
      dataset_id = rep(i,length(trans)),
      imp_id = rep(j,length(trans)),
      trans = trans,
      rate = rate,
      rate_se = rate_se,
      shape = shape,
      shape_se = shape_se,
      beta_educ = beta_educ,
      beta_educ_se = beta_educ_se,
      beta_sex = beta_sex,
      beta_sex_se = beta_sex_se
    )
  }else {
    est_obj <- data.frame(
    model = rep(model_name, length(trans)),
    dataset_id = rep(i,length(trans)),
    trans = trans,
    rate = rate,
    rate_se = rate_se,
    shape = shape,
    shape_se = shape_se,
    beta_educ = beta_educ,
    beta_educ_se = beta_educ_se,
    beta_sex = beta_sex,
    beta_sex_se = beta_sex_se
  )
  }
  
  return(est_obj)
}


rename_trans<- function(names){
  
  new_names <- gsub("Transition (\\d+) to (\\d+)", "\\1->\\2", names)
  return(new_names)
}
