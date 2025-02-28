library(magrittr)
library(tidyverse)
library(gtsummary)
library(ggplot2)
library(ggprism)
library(ggpubr)
library(poLCA)
library(sBIC)
library(dplyr) 
library(future.apply)
library(parallel)
library(tictoc)
library(purrr)

source("Functions//flexsurv.R")
scenario <- "B"
sim_obj_name <- load(paste0("Simulation/data for simulation/scenario_",scenario,".RData"))
sim_obj <- get(sim_obj_name)
colnames(sim_obj$pattern_obj$obj$y)<- tolower(colnames(sim_obj$pattern_obj$obj$y))

nsim <- 10000
N<- 100
load(paste0("Simulation/Simulation Outputs/pop_study_",scenario,"_",N,"_",nsim,".RData"))

X <- pop %>%  dplyr::select(any_of(colnames(sim_obj$pattern_obj$obj$y))) %>% dplyr::mutate_all(function(x)x+1)
if (scenario %in% c("A", "Av2")){
  post <- poLCA.posterior(sim_obj$pattern_obj$obj,y=X)
} else if (scenario %in% c("B")){
  post <- poLCA.posterior(sim_obj$pattern_obj$obj,y=X,x=as.matrix(cbind(pop$age), ncol=1))
}
post

## auxiliary functions

create_unique_folder <- function(prefix = "results") {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  folder_name <- paste0(prefix, "_", timestamp,"_imp")
  return(folder_name)
}
run_multistate <- function(pop,result_folder){
  
  pop_ms <- pop %>%
    mutate(death = ifelse(Age_exit < Age_death, 0, 1))
  pop_ms_death_1 <- pop_ms %>%
    dplyr::filter(death == 1) %>%
    group_by(dataset_id, patient_id) %>%
    group_split() %>%
    map_df(~ bind_rows(
      .x,
      slice_tail(.x, n = 1) %>% mutate(MP = dim(sim_obj$tmat)[1],MP_sim=dim(sim_obj$tmat)[1], age = Age_death)
    )) %>%
    ungroup()
  
  pop_ms_death_0 <- pop_ms %>%
    dplyr::filter(death == 0)
  
  pop_ms <- bind_rows(pop_ms_death_0, pop_ms_death_1)
  
  pop_ms <- pop_ms %>%
    group_by(dataset_id,patient_id) %>%
    filter(n() > 1) %>%
    ungroup()
  
  if (scenario=="B"){
    pop_ms <- pop_ms %>%
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
  apply_flexsurv_imp(pop_ms, sim_obj, result_folder)
}
run_multistate_wrapper <- function(dataset) {
  run_multistate(dataset, result_folder)
}

######## imputation ############
set.seed(123)

n_imputations <- 100
n_batch <- 10
el_per_batch <- 10 # n_batch*el_per_batch should be equal to n_imputations
tic("Imputation")
MP_imputed <- apply(post, 1, function(p) {
  sample(1:ncol(post), size = n_imputations, replace = TRUE, prob = p)
})
MP_imputed <- t(MP_imputed)


MP_prev <- data.frame(p1=rowSums(MP_imputed==1)/ncol(MP_imputed),p2=rowSums(MP_imputed==2)/ncol(MP_imputed))
diff <- data.frame(d1 = abs(MP_prev$p1-post[,1]), d2=abs(MP_prev$p2-post[,2]))

result_folder <-  file.path("results",create_unique_folder(paste0("results")))
dir.create(result_folder, recursive = TRUE)

column_names <- c("dataset_id", "patient_id"  ,  "Age_death"  ,   "Age_entry"   ,  "educ_el"     ,  "dm_sex"   ,     "time_in_study" ,"base_MP"   ,   
                  "ndis"       ,   "Age_exit"  ,    "MP_base"  ,     "age"      ,     "visit_number" , "age_group"   ,  "lag_age"  ,     "MP_sim"   ,    
                  "MP"         ,   "imp_nr")

load("Data/disease_names.RData")
pop <- pop[, setdiff(colnames(pop), disease_names)]

cores <- min(detectCores()-1, 24)
plan(multisession, workers=cores)  # Or `plan(multiprocess)` for cross-platform compatibility

for (j in 1:n_batch){
  imputed_datasets <- do.call(rbind, lapply(1:el_per_batch, function(i) {
    index <- (j - 1) * el_per_batch + i
    if (index > n_imputations) {
      print("Index exceeds n_imputations!")
      return(NULL)
    }
    temp_pop <- pop
    temp_pop$MP <- MP_imputed[, index]
    temp_pop$imp_nr <- index
    return(temp_pop)
  }))
  if (is.null(imputed_datasets) || nrow(imputed_datasets) == 0) {
    message("Skipping batch ", j, " as no valid data was generated.")
    next  
  }
  #imputed_datasets <- do.call(rbind, imputed_datasets[!sapply(imputed_datasets, is.null)])
  datasets <- split(imputed_datasets, interaction(imputed_datasets$dataset_id, imputed_datasets$imp_nr))

  # apply flexsurv to imputed datasets
  future_lapply(datasets, FUN = run_multistate_wrapper)
  
}
toc()
#####


