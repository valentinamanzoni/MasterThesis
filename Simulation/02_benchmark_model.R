
library(tidyverse)
library(magrittr)
library(dplyr)
library(poLCA)
library(parallel)
library(future.apply)
library(tictoc)
############# EXACT TIMES ########
#################################
# implement dataset with exact times of transition (from one cluster to another)

N <- 100
scenario <- "B"
# how many individuals?
nsim <- 10000

sim_obj_name <- load(paste0("Simulation/data for simulation/scenario_",scenario,".RData"))
sim_obj <- get(sim_obj_name)
colnames(sim_obj$pattern_obj$obj$y)<- tolower(colnames(sim_obj$pattern_obj$obj$y))

# import the simulated data with exact time of transitions and exact MP
load(paste0("Simulation/Simulation Outputs/sim_MP_",scenario,"_",N,"_",nsim,".RData"))

# add Age_exit
sim_MPB <- sim_MPB %>% 
  group_by(dataset_id,patient_id)%>%
  mutate(Age_exit=min(Age_death, Age_entry + time_in_study)) %>%
  ungroup()

sim_MPB <- sim_MPB %>% 
  relocate(Age_exit, .after = Age)


sim_MPB <- sim_MPB %>%
  mutate(death = ifelse(Age_exit < Age_death, 0, 1))
pop_ms_death_1 <- sim_MPB %>%
  mutate(MP = if_else(death == 1 & Age == Age_death, dim(sim_obj$tmat)[1], MP))

# rename MP of the simulation (ground truth) to MP_sim 
pop_ms_death_1 <- pop_ms_death_1%>%
  rename(MP_sim=MP, age=Age)

#
source("Functions//Apply_LCA.R")
pop<- apply_LCA(pop_ms_death_1,sim_obj, scenario)
pop <- pop %>%
  mutate(MP = if_else(death == 1 & age == Age_death, dim(sim_obj$tmat)[1], MP))

source("Functions//flexsurv.R")
# auxiliary functions:
create_unique_folder <- function(prefix = "results") {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  folder_name <- paste0(prefix, "_", timestamp,"_imp")
  return(folder_name)
}

## fit the benchmark model on exact time of transitions and exact states
run_multistate_ex <- function(pop,result_folder){
  
  if (scenario=="B"){
    pop <- pop %>%
      group_by(patient_id) %>%
      mutate(flag = cumsum(MP == 2),
             # Adjust MP based on the flag
             MP = case_when(
               flag > 0 & MP == 1 ~ 2, # If flag is active and MP is 1, change it to 2
               TRUE ~ MP           # Otherwise, keep MP as it is
             )
      ) %>%
      ungroup() %>%
      dplyr::select(-flag) 
  }
  apply_flexsurv(pop, sim_obj, result_folder)
  apply_flexsurv_base(pop, sim_obj, result_folder)
  gc()
  
}
run_multistate_wrapper_ex <- function(dataset) {
  run_multistate_ex(dataset, result_folder)
}

cores <- detectCores()-5 # modify number of cores
plan(multisession, workers=cores)  # Or `plan(multiprocess)` for cross-platform compatibility

datasets <- split(pop, pop$dataset_id)

# create result folder using timestamp
result_folder <-  file.path("results",create_unique_folder(paste0("results")))
dir.create(result_folder, recursive = TRUE)


tic("Multistate in parallel")
future_lapply(datasets, FUN = run_multistate_wrapper_ex)
toc()



########################################
######## Estimated Parameters###########
#######################################
# extract the parameters estimated by the bechmark model
library(ggplot2)
source("Functions//extract_performance.R")


# folder where estimates will be saved
result_folder <-  "results/results_estimates"
# folder from which the fitted models are imported
results_path <- paste0("results/results_20250120_121006_imp/")


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
models<-c("flexsurv","flexsurv_base")
load_model <- FALSE
if (load_model){
  pt <-load(paste0("results/results_estimates/model_est_ET_",scenario,"_",nsim,".RData"))
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
model_est <-list(
  param_df = param_df,
  time_df = conv_time
)

if (!load_model){
  save(model_est, file =paste0(result_folder,"/model_est_ET_",scenario,"_", nsim,".RData"))
}

