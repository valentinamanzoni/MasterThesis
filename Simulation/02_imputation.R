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

nsim <- 3000
N<- 100
load(paste0("Simulation/Simulation Outputs/pop_study_",scenario,"_",N,"_",nsim,".RData"))

X <- pop %>%  dplyr::select(any_of(colnames(sim_obj$pattern_obj$obj$y))) %>% dplyr::mutate_all(function(x)x+1)
if (scenario %in% c("A", "Av2")){
  post <- poLCA.posterior(sim_obj$pattern_obj$obj,y=X)
} else if (scenario %in% c("B")){
  post <- poLCA.posterior(sim_obj$pattern_obj$obj,y=X,x=as.matrix(cbind(pop$age), ncol=1))
}
post

######## imputation ############
set.seed(123)
n_imputations <- 100
MP_imputed <- apply(post, 1, function(p) {
  sample(1:ncol(post), size = n_imputations, replace = TRUE, prob = p)
})
MP_imputed <- t(MP_imputed)

imputed_datasets <- do.call(rbind, lapply(1:n_imputations, function(i) {
  temp_pop <- pop
  temp_pop$MP <- MP_imputed[, i]
  temp_pop$imp_nr <- i  
  return(temp_pop)
}))
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
# split across imputed datasets
run_multistate_wrapper <- function(dataset) {
  run_multistate(dataset, result_folder)
}

######
result_folder <-  file.path("results",create_unique_folder(paste0("results")))
dir.create(result_folder, recursive = TRUE)

output_file <- file.path(result_folder, "imputed_datasets.csv")

column_names <- c("dataset_id", "patient_id"  ,  "Age_death"  ,   "Age_entry"   ,  "educ_el"     ,  "dm_sex"   ,     "time_in_study" ,"base_MP"   ,   
                  "ndis"       ,   "Age_exit"  ,    "MP_base"  ,     "age"      ,     "visit_number" , "age_group"   ,  "lag_age"  ,     "MP_sim"   ,    
                   "MP"         ,   "imp_nr")

empty_df <- data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(empty_df) <- column_names


write.table(empty_df, file = output_file, col.names = TRUE, row.names = FALSE, sep = ",")


load("Data/disease_names.RData")
pop <- pop[, setdiff(colnames(pop), disease_names)]
head(pop)
for (i in 1:n_imputations) {
  temp_pop <- pop
  temp_pop$MP <- MP_imputed[, i]
  temp_pop$imp_nr <- i
  
  write.table(temp_pop, file = output_file, append = TRUE, col.names = FALSE, row.names = FALSE, sep = ",")
  print(i)
}

#disease_names <- colnames(pop)[17:ncol(pop)]
#save(disease_names, file = "Data/disease_names.RData")

# load Imputed_datasets
imputed_datasets <- read.csv(paste0(result_folder,"/imputed_datasets.csv"), sep = ",",stringsAsFactors = FALSE, header = FALSE)
head(imputed_datasets)
colnames(imputed_datasets)<- column_names
imputed_datasets <- imputed_datasets%>%
  mutate(dataset_id=as.integer(dataset_id),
         patient_id=as.integer(patient_id),
         educ_el=as.integer(educ_el),
         dm_sex=as.integer(dm_sex),
         ndis=as.integer(ndis),
         MP= as.numeric(MP),
         MP_sim= as.numeric(MP_sim),
         age= as.numeric(age),
         age= as.numeric(age),
         Age_death= as.numeric(Age_death),
         Age_entry= as.numeric(Age_entry),
         Age_exit= as.numeric(Age_exit)
         )

cores <- min(detectCores(), 24)
plan(multisession, workers=cores)  # Or `plan(multiprocess)` for cross-platform compatibility


datasets <- split(imputed_datasets, interaction(imputed_datasets$dataset_id, imputed_datasets$imp_nr))


# apply flexsurv to imputed datasets
tic("Imputation")
future_lapply(datasets, FUN = run_multistate_wrapper)
toc()

#####

# drop diseases to reduce memory usage:
# diseases <-colnames(sim_obj$pattern_obj$obj$y)
# imputed_datasets <- imputed_datasets[, setdiff(colnames(imputed_datasets), diseases)]
# 
# cores <- detectCores()
# plan(multisession, workers=cores)  # Or `plan(multiprocess)` for cross-platform compatibility
# 
# 
# datasets <- split(imputed_datasets, interaction(imputed_datasets$dataset_id, imputed_datasets$imp_nr))
# 
# 
# # apply flexsurv to imputed datasets
# tic("Imputation")
# future_lapply(datasets, FUN = run_multistate_wrapper)
# toc()
# 
# 
