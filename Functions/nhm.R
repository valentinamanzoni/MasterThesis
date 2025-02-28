library(nhm)
library(msm)
library(tictoc)

# result <- matrix(0, nrow = nrow(misc), ncol = ncol(misc))
# result[which(misc!= 0)] <- seq_along(which(misc != 0))

apply_nhm <- function(pop_ms,misc, sim_obj, result_folder){
  if (scenario %in% c("A", "Av2")){
    #pop_ms$state <- as.numeric(as.character(pop_ms$MP))
    pop_ms$state <- pop_ms$MP
    pop_ms <- as.data.frame(pop_ms)
    q_matrix <- sim_obj$tmat
    colnames(q_matrix)<-c(1:dim(q_matrix)[1])
    rownames(q_matrix)<-c(1:dim(q_matrix)[1])
    q_matrix[is.na(q_matrix)] <- 0
    # same structure for q_mat and nonh
    nonh <- q_matrix
    tic("nhm gompertz")
    model_obj_gomp<- model.nhm(state ~ age, subject=patient_id, type='gompertz', data=pop_ms, trans= q_matrix,
                           nonh= nonh,
                           #initp = c(1,2,2,2,2,0),
                           #centre_time=72,   
                           death=T, death.states = 6)

    model_0 <- nhm(model_obj_gomp,
                   #initial= nhm_init,
                   gen_inits = TRUE, #not converging
                   control=nhm.control(splits = c(60,61,63,65,75,76,85,95,99,102,103,105,109), verbose=TRUE)
                   )
    t<-toc()
    
    model_obj_0 <- list(
      model = model_0,
      time = t$toc - t$tic
    )
    save(model_obj_0, file =paste0(result_folder,"/nhm_base_scenario_",scenario,"_", nsim,"_",unique(pop_ms$dataset_id),".RData"))
    
    
    non_diagonal <- misc
    diag(non_diagonal) <- 0
    misc_param <- non_diagonal[non_diagonal != 0]
    
    # initialize parameters for misclassified model
    nhm_init <- c(model_0$par, misc_param)

    misc_shape <- matrix(0, nrow = nrow(non_diagonal), ncol = ncol(non_diagonal))
    misc_shape[which(non_diagonal!= 0)] <- seq_along(which(non_diagonal != 0)) # numera da 1 length misc_param gli stati non nulli della matrice non_diag
      
    tic("nhm misc gompertz")
    model_obj_gomp_misc<- model.nhm(state ~ age, subject=patient_id, type='gompertz', data=pop_ms, trans= q_nhm,
                                    nonh= nonh,
                                    emat = misc_shape,
                                    #initp = c(1,2,2,2,2,0),
                                    # initp_value
                                    #centre_time=72,   
                                    death=T, death.states = 3)
    
    # misc parameters need to be passed to initial and then fixed in fixedpar
    n1<- length(model_0$par) + 1
    n2<- length(nhm_init)
    model_misc <- nhm(model_obj_gomp_misc,
                      initial= nhm_init,
                      #gen_inits = TRUE,
                      fixedpar= c(n1:n2),
                      control=nhm.control(splits = c(60,61,63,65,75,76,85,95,99,102,103,105,109)
                                          , verbose=TRUE
                      )
    )
    t<-toc()
    model_obj_m <- list(
      model = model_misc,
      time = t$toc - t$tic
    )
    save(model_obj_m, file =paste0(result_folder,"/nhm_misc_scenario_",scenario,"_", nsim,"_",unique(pop_ms$dataset_id),".RData"))
    
    
    
  } else if (scenario=="B"){
    model_0 <- NULL
    model_misc <- NULL
    pop_ms$state <- pop_ms$MP
    pop_ms <- as.data.frame(pop_ms)
    q_nhm <- rbind(c(0,1,2),
                   c(0,0,3),
                   c(0,0,0))
    colnames(q_nhm)<-c(1:dim(q_nhm)[1])
    rownames(q_nhm)<-c(1:dim(q_nhm)[1])
    # same structure for q_mat and nonh
    nonh <- q_nhm
    covm <- list(
      educ_el = rbind(c(0,1,2), c(0,0,3), c(0,0,0)),
      dm_sex = rbind(c(0,4,5), c(0,0,6), c(0,0,0)))
    tic("nhm gompertz with cov")
    model_obj_gomp<- model.nhm(state ~ age, subject=patient_id, type='gompertz', data=pop_ms, trans= q_nhm,
                               nonh= nonh,
                               covariates= c("educ_el", "dm_sex"),
                               covm = covm,
                               #initp = c(1,2,2,2,2,0),
                               # initp_value
                               #centre_time=72,   
                               death=T, death.states = 3)
    tryCatch({
    model_0 <- nhm(model_obj_gomp,
                   #initial= nhm_init,
                   gen_inits = TRUE, #not converging
                   control=nhm.control(splits = c(60,61,63,65,75,76,85,95,99,102,103,105,109),
                                       ncores = 4
                                       #, verbose=TRUE
                   )
    )
    }, error = function(e) {
      print(paste("Error during model fitting:", e$message))
    })
    
    t<-toc()
    
    if (!is.null(model_0)){
      model_obj_0 <- list(
        model = model_0,
        time = t$toc - t$tic
      )
      save(model_obj_0, file =paste0(result_folder,"/nhm_base_scenario_",scenario,"_", nsim,"_",unique(pop_ms$dataset_id),".RData"))
    } else {
      cat("Model nhm_base is NULL for dataset_id ",unique(pop_ms$dataset_id),". Model was not saved", "\n")
     }
    
    non_diagonal <- misc
    diag(non_diagonal) <- 0
    misc_param <- non_diagonal[non_diagonal != 0]
    
    #nhm_init <- c(model_0$par, misc_param)
    nhm_init <- c(model_0$par, log(misc_param))
    
    misc_shape <- rbind(c(0,2,0),
                        c(1,0,0),
                        c(0,0,0))
    tic("nhm misc gompertz with cov")
    model_obj_gomp_misc<- model.nhm(state ~ age, subject=patient_id, type='gompertz', data=pop_ms, trans= q_nhm,
                                    nonh= nonh,
                                    emat = misc_shape,
                                    covariates= c("educ_el", "dm_sex"),
                                    covm = covm,
                                    #firstobs="exact",
                                    #initp = c(1,2,2,2,2,0),
                                    #initp_value = c(0.6,0.4,0),
                                    #centre_time=72,   
                                    death=T, death.states = 3)
    
    
    # misc parameters need to be passed to initial and then fixed in fixedpar
    n1<- length(model_0$par) + 1
    n2<- length(nhm_init)
    split_points <- c(60,61,63,65,75,76,85,95,99,102,103,105,106,109)
    # split_points <- find_splits(pop_ms$age)
    #model_misc <- recursive_nhm(model_obj_gomp_misc,nhm_init,c(n1:n2), split_points, 20)
    
    tryCatch({
      model_misc <- nhm(
        model_obj_gomp_misc,
        initial = nhm_init,
        fixedpar = c(n1:n2),
        control = nhm.control(
          splits = split_points,
          ncores = 4, 
          rtol = 1e-8, 
          atol = 1e-8
        )
      )
    }, error = function(e) {
      print(paste("Error during model fitting:", e$message))
    })
    t<-toc()
    if (!is.null(model_misc)) {
      model_obj_m <- list(
        model = model_misc,
        time = t$toc - t$tic
      )
      save(model_obj_m, file =paste0(result_folder,"/nhm_misc_scenario_",scenario,"_", nsim,"_",unique(pop_ms$dataset_id),".RData"))
      
    } else {
      cat("Model nhm_misc is NULL for dataset_id ",unique(pop_ms$dataset_id),". Model was not saved")
    
    }
    
  }
  else
    stop("Insert a valid design type")
  print("nhm executed")
}


find_splits <- function(age) {
  quantiles <- seq(min(age),max(age), length.out=100)
  #quantiles <- quantile(age, probs = seq(0, 1, 0.005))
  return(quantiles[-c(1, length(quantiles))])
}

recursive_nhm <- function(model_obj_gomp_misc,nhm_init,fixed_param, split_points,c){
  tryCatch({
    cat("Updated split points: ",split_points)
    cat("nhm_init",nhm_init)
    model_misc <- nhm(
      model_obj_gomp_misc,
      initial = nhm_init,
      fixedpar = c(n1:n2),
      control = nhm.control(
        splits = split_points,
        ncores = 4, 
        rtol = 1e-6, 
        atol = 1e-6
      )
    )
    return(model_misc) # stop if successful
  }, error = function(e) {
    print(paste("Error during model fitting:", e$message))
    
    if (grepl("Matrix singular when inverting transition probability matrix", e$message)) {
      problematic_points <- as.numeric(unlist(regmatches(e$message, gregexpr("[0-9]+\\.[0-9]+", e$message))))
      if (length(problematic_points) == 2) {
        # Add a new split point between the problematic points
        new_split <- mean(problematic_points)
        split_points <- sort(unique(c(split_points, new_split)))
        print(split_points)
        print(paste("Added new split at:", new_split))
        if (c>=0)
          return(recursive_nhm(model_obj_gomp_misc,nhm_init,fixed_param, split_points,c-1))
        else 
          return(NULL)
      } else {
        stop("Unable to parse problematic points from error message.")
      }
    } else {
      print(paste("An unexpected error occurred:", e$message))
      return(NULL)
    }
  })
}


# 
# error <- F
# 
# 
# tryCatch({
#     object_nhm <- model.nhm(state ~ age,
#                             subject = patient_id,
#                             data = temp,
#                             trans = tmat_1,
#                             nonh = tmat_1,
#                             type = "gompertz",
#                             covariates = c("cov1", "cov2", "cov3"),
#                             covm = list(cov1= tmat_1, cov2=tmat_2, cov3=tmat_3),
#                             censor.states = c(1,2),
#                             censor = 99,
#                             death = T,
#                             death.states = c(3))
#     
#     
#     model_nhm <- nhm(object_nhm,
#                      gen_inits = T,
#                      #initial = initial_guess,
#                      score_test = FALSE,
#                      control = nhm.control(ncores = cores_nhm, obsinfo = FALSE, coarsen = T, coarsen.vars = c(1), coarsen.lv = 5, splits = split_points, rtol=1e-4, atol=1e-4))
#   },
# error = function(e) {
#     print(paste("Error during model fitting:", e$message))
#     error <<- TRUE
#   })

