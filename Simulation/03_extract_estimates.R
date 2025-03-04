library(flexsurv)
library(flexsurvcure)
library(ggplot2)
library(dplyr)
library(survminer)
library(tidyverse) 


###########################################################################################
#################### Extract estimates and performances from fitted model #################
###########################################################################################

source("Functions//extract_performance.R")

## load desired models
nsim<-10000
N <-100 # total number of datasets

# study type: "ps","rv", "ET"
st <- "rv"

# folder from where
results_path <- paste0("results/results_20250120_160025_rv","/")
scenario <- "B"
# folder where estimates will be saved
result_folder <-  "results/results_estimates"

if (!dir.exists(result_folder)) {
  dir.create(result_folder)
  message("Folder 'results_estimates' created.")
} else {
  message("Folder 'result_estimate' already exists.")
}
models<-c("flexsurv","flexsurv_base", "nhm_base","nhm_misc","model_age", "model_age_misc")
#models<-c("flexsurv","flexsurv_base")

## load true parameters
gt <-load(paste0("Simulation/data for simulation/scenario_",scenario,".RData"))
true_param_obj <- get(gt)
true_param_obj$trans_mod
ntrans <- length(true_param_obj$trans_mod)
true_param<- data.frame(trans=character(),
                        rate=numeric(),
                        shape=numeric(),
                        beta_educ=numeric(),
                        beta_sex = numeric())
trans_label <- c("1->2","1->3","2->3")
for (i in 1:ntrans){
  curr <- true_param_obj$trans_mod[[i]]
  true_param <-rbind(true_param, data.frame(trans=trans_label[i], # DA CAMBIARE
                                        rate= curr$res.t["rate", "est"],
                                        shape= curr$res.t["shape", "est"],
                                        beta_educ= curr$res.t["educ_el", "est"],
                                        beta_sex= curr$res.t["dm_sexFemales", "est"]
                                        ))
}
true_param
#save(true_param, file =paste0(result_folder,"/true_param.RData"))


########## to just add new models ########
load_model <- FALSE
if (load_model){
    pt <-load(paste0("results/results_estimates/model_est_",scenario,"_",nsim,".RData"))
    param <- get(pt)
    param_df <- param$param_df
    conv_time <- param$time_df
} else {
  param_df <-data.frame(dataset_id = integer(), 
                        trans=character(),
                        rate=numeric(),
                        shape=numeric(),
                        rate_se=numeric(),
                        shape_se=numeric(),
                        beta_educ = numeric(),
                        beta_educ_se = numeric(),
                        beta_sex = numeric(),
                        beta_sex_se = numeric()
                        )
  
  # extract rate, shape, their se and the transition name and store it in param for each dataset
  conv_time <- data.frame(dataset_id = integer(), 
                          model_name = character(),
                          time = numeric()
  )
  # param_df <- param_df %>% filter(model!="model_age_misc")
  # conv_time <- conv_time%>%filter(model_name!="model_age_misc")
  for(model_name in models){
    print(model_name)
    for(i in 1:N){
      model_file <- paste0(results_path, model_name, "_scenario_", scenario, "_", nsim, "_", i, ".RData")
      model_n <- tryCatch({
        load(model_file)
      }, error = function(e) {
        message(paste("File not found or cannot load:", model_file))
        return(NULL) 
      })
      
      if (is.null(model_n)) {
        next
      }
      model_n <-load(paste0(results_path,model_name,"_scenario_",scenario,"_",nsim,"_",i,".RData"))
      model <- get(model_n)
      conv_time <- rbind(conv_time, data.frame(dataset_id = i, model_name = model_name, time=model$time))
      #print(param_df)
      param_df <- rbind(param_df,get_estimates(model, model_name,i))
  }
  }
}
param_df
conv_time
param_df$model <- dplyr::recode(param_df$model,
                                "flexsurv_base" = "flexsurv_true",
                                "model_age_misc" = "msm_misc",
                                "model_age" = "msm_base"
                                )
conv_time$model_name <- dplyr::recode(conv_time$model_name,
                                 "flexsurv_base" = "flexsurv_true",
                                 "model_age_misc" = "msm_misc",
                                 "model_age" = "msm_base"
)
###### Imputation param
imp_param <-FALSE
n_imp <- 20
results_path <- paste0("results/","results_20241224_171609","/")
model_name <- "flexsurv_imp"
if (imp_param){
  param_df_imp <-data.frame(dataset_id = integer(), 
                        imp_id = integer(), 
                        trans=character(),
                        rate=numeric(),
                        shape=numeric(),
                        rate_se=numeric(),
                        shape_se=numeric(),
                        beta_educ = numeric(),
                        beta_educ_se = numeric(),
                        beta_sex = numeric(),
                        beta_sex_se = numeric()
  )
  
    conv_time_imp <- data.frame(dataset_id = integer(), 
                          imp_id = integer(), 
                          model_name = character(),
                          time = numeric())
 
  for(j in 1:n_imp){
    for(i in 1:N){
      model_file <- paste0(results_path, model_name, "_scenario_", scenario, "_", nsim, "_", i,"_",j, ".RData")
      model_n <- tryCatch({
        load(model_file)
      }, error = function(e) {
        message(paste("File not found or cannot load:", model_file))
        return(NULL) 
      })
      
      if (is.null(model_n)) {
        next
      }
      #model_n <-load(paste0(results_path,model_name,"_scenario_",scenario,"_",nsim,"_",i,".RData"))
      model <- get(model_n)
      conv_time_imp <- rbind(conv_time_imp, data.frame(dataset_id = i, imp_id = j, model_name = model_name, time=model$time))
      #print(param_df)
      param_df_imp <- rbind(param_df_imp,get_estimates(model, model_name,i,j))
    }
  }
  
    conv_time2 <- conv_time_imp%>%
      group_by(dataset_id) %>%
      mutate(model_name=model_name,
        time= sum(time)) %>%
      dplyr::select(-imp_id)
    #view(conv_time2)
    param_df3 <- param_df_imp %>%
      group_by(dataset_id, trans) %>%
      summarise(
        rate_se_imp = sqrt(mean(rate_se^2, na.rm = TRUE) + (1 + 1/n_imp) * var(rate, na.rm = TRUE)),
        shape_se_imp = sqrt(mean(shape_se^2, na.rm = TRUE) + (1 + 1/n_imp) * var(shape, na.rm = TRUE)),
        beta_educ_se_imp = sqrt(mean(beta_educ_se^2, na.rm = TRUE) + (1 + 1/n_imp) * var(beta_educ, na.rm = TRUE)),
        beta_sex_se_imp = sqrt(mean(beta_sex_se^2, na.rm = TRUE) + (1 + 1/n_imp) * var(beta_sex, na.rm = TRUE)),
        
        rate = mean(rate, na.rm = TRUE),
        shape = mean(shape, na.rm = TRUE),
        beta_educ = mean(beta_educ, na.rm = TRUE),
        beta_sex = mean(beta_sex, na.rm = TRUE),
        
        .groups = "drop" 
             )
    
    param_df_new <-data.frame(model=rep("flexsurv_imp",length(param_df3$dataset_id)),
                              dataset_id = param_df3$dataset_id, 
                              trans=param_df3$trans,
                              rate=param_df3$rate,
                              rate_se=param_df3$rate_se_imp,
                              shape=param_df3$shape,
                              shape_se=param_df3$shape_se_imp,
                              beta_educ = param_df3$beta_educ,
                              beta_educ_se = param_df3$beta_educ_se_imp,
                              beta_sex = param_df3$beta_sex,
                              beta_sex_se = param_df3$beta_sex_se_imp
    )
    param_df <- param_df%>% filter(model!="flexsurv_imp")
    param_df <- rbind(param_df,param_df_new)
    conv_time <- rbind(conv_time,conv_time2)
    


    }


model_est <-list(
  param_df = param_df,
  time_df = conv_time
)

if (!load_model){
  save(model_est, file =paste0(result_folder,"/model_est_",st,"_",scenario,"_", nsim,".RData"))
}

######

######
######################################
######## Estimates Plot ###########
########################################
#use rate and shape, covariates parametes and their se,to produce the plots.
source("Functions//plot_estimates.R")
boxplot_param(param_df,true_param,nsim)
dotplot_param(param_df,true_param,nsim)
#plot_param_hist(param_df,true_param,nsim)
plot_list <-plot_param_hist_2(param_df,true_param,nsim)
plot_list$rate
plot_list$shape
plot_list$beta_educ
plot_list$beta_sex
########################################
######## Performance measures ###########
########################################
# compute performance measures: bias, empirical standard error, MSE, Coverage, bias-eliminated coverage
source("Functions//performance_measures.R")
pm_df <- compute_pm(param_df,true_param) 
#pm_df <- compute_pm2(param_df,true_param)
#view(pm_df)
# plot relative biases
#plot_rel_bias(pm_df)
dotplot_rel_bias(param_df, true_param,nsim)
if (!load_model){
  save(pm_df, file =paste0(result_folder,"/pm_df_",scenario,"_", nsim,".RData"))
}
#########################################
######## Convergence Analysis ###########
#########################################
# estimate average convergence time 
boxplot(time~ model_name,conv_time)


##########################################
######### Other Plots ####################
#######################################
