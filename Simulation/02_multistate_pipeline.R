
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

# change the following to load correct scenario and simulated data
# load simulated data
scenario <- "B"
study_type <- "ps" # random visit rv or pop study ps
sim_obj_name <- load(paste0("Simulation/data for simulation/scenario_",scenario,".RData"))
sim_obj <- get(sim_obj_name)
colnames(sim_obj$pattern_obj$obj$y)<- tolower(colnames(sim_obj$pattern_obj$obj$y))

nsim <- 10000
N<- 100

if (study_type=="ps"){
  data<-load(paste0("Simulation/Simulation Outputs/pop_study_",scenario,"_",N,"_",nsim,".RData"))
  pop <- get(data)
} else if (study_type=="rv"){
  data <- load(paste0("Simulation/Simulation Outputs/rv_pop_study_",scenario,"_",N,"_",nsim,".RData"))
  pop <- get(data)
} else{
  cat("Insert valid study type")
}

#pop <- pop %>% filter(dataset_id %in% c(91))


######### apply LCA ##########
# load LCA parameter & use them to assign MP to the data
source("Functions//Apply_LCA.R")

pop<- apply_LCA(pop,sim_obj, scenario)
head(pop)
############## Multistate ################
#######################################
source("Functions//misc_matrix.R")
# add check of q_matrix at 1sr, 2nd and 3rd follow up per dataset
datasets <- split(pop, pop$dataset_id)
m_matrices <- lapply(datasets, calculate_m_matrix, index = 1)
m_matrices_2 <- lapply(datasets, calculate_m_matrix, index = 2)
m_matrices_3 <- lapply(datasets, calculate_m_matrix, index = 3)

m_mat_approx <-lapply(datasets, calc_approx_m_matrix, index = 1)
  
n <- dim(sim_obj$tmat)[1]
m_matrix_array <- array(unlist(m_matrices), dim = c(n, n, length(m_matrices)))
m_matrix_array_2 <- array(unlist(m_matrices_2), dim = c(n, n, length(m_matrices_2)))
m_matrix_array_3 <- array(unlist(m_matrices_3), dim = c(n, n, length(m_matrices_3)))
m_approx_array  <- array(unlist(m_mat_approx), dim = c(n, n, length(m_mat_approx)))

# mean and sample variance of m_matrix
print_mean_m_mat <- function(m_matrix_array){
  mean_matrix <- apply(m_matrix_array, c(1, 2), mean)
  rownames(mean_matrix) <- rownames(sim_obj$tmat)
  colnames(mean_matrix) <- rownames(sim_obj$tmat)
  variance_matrix <- apply(m_matrix_array, c(1, 2), var)
  rownames(variance_matrix) <- rownames(sim_obj$tmat)
  colnames(variance_matrix) <- rownames(sim_obj$tmat)
  print("mean misc matrix")
  print(mean_matrix)
  print("variance misc matrix")
  print(variance_matrix)
}

print_mean_m_mat(m_matrix_array)
print_mean_m_mat(m_matrix_array_2)
print_mean_m_mat(m_matrix_array_3)
print_mean_m_mat(m_approx_array)

#############################
source("Functions//msm.R")
source("Functions//nhm.R")
source("Functions//flexsurv.R")



















# auxiliary functions:
create_unique_folder <- function(prefix = "results") {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  folder_name <- paste0(prefix, "_", timestamp,"_",study_type)
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
  
  misc <- calculate_m_matrix(pop,1)
  apply_flexsurv(pop_ms, sim_obj, result_folder)
  apply_msm(pop_ms,misc, sim_obj, result_folder)
  apply_flexsurv_base(pop_ms, sim_obj, result_folder)
  gc()
  apply_nhm(pop_ms,misc, sim_obj,result_folder)
  gc()
  
}
run_multistate_wrapper <- function(dataset) {
  run_multistate(dataset, result_folder)
}


########### split by dataset_id #########
# run in parallel
# group pop by dataset_id
cores <- min(detectCores() -5, 24 ) #modify number of cores
plan(multisession, workers=cores)  # Or `plan(multiprocess)` for cross-platform compatibility

datasets <- split(pop, pop$dataset_id)

# create result folder using timestamp
result_folder <-  file.path("results",create_unique_folder(paste0("results")))
dir.create(result_folder, recursive = TRUE)


# call function to apply nhm, msm, and flexsurv in parallel
tic("Multistate in parallel")
future_lapply(datasets, FUN = run_multistate_wrapper)
toc()

