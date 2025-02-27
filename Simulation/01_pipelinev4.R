library(tidyverse)
library(magrittr)
library(flexsurv)
library(mstate)
library(hesim)
library(data.table)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(splines)
library(stats)
library(rstpm2)
library(truncdist)
library(flexsurv)
source("Simulation//fun//fun_simulate_MM_trajv4.R")
source("Simulation//fun//Study_design_fun_2.R")

# set you folder path
path_sim <- "Simulation"

###########################################################
######   SCENARIO AV2 #######
###########################################################

# 1) How many datasets?
N <- 100

# 6) how many individuals?
nsim <- 100000

scenario <- "Av2"
load(paste0("Simulation/data for simulation/scenario_",scenario,".RData"))

sim_test_list <- sim_mm_traj(nsim,N,scenario_obj = scenario_obj_Av2,seed=1258,scenario = "Av2")
sim_test <- sim_test_list[[1]]
sim_MP <- sim_test_list[[2]]

# save simulation results 
output_folder <-paste0(path_sim, "/", "Simulation Outputs")

if (!dir.exists(output_folder)) {
  dir.create(output_folder)
  message("Folder 'Simulation Outputs' created.")
} else {
  message("Folder 'Simulation Outputs' already exists.")
}

save(sim_test, file = paste0(output_folder, "/", "sim_data_",scenario,"_",as.character(N),"_",as.character(nsim),".RData"))
save(sim_MP, file = paste0(output_folder, "/", "sim_MP_",scenario,"_",as.character(N),"_",as.character(nsim),".RData"))
message("Dataset saved in 'Simulation Outputs' folder.")


design_type <- "popstudy"

####### CALL TO STUDY DESIGN FUNCTION ####
if (design_type =="popstudy"){
  age_tr <- 78
  obs_t_int <- c(6,3)
  pop <- popbased_study(sim_test,sim_MP, underreporting = FALSE,
                             obs_t_int = obs_t_int, age_tr = age_tr)
  
  save(pop, file = paste0(output_folder, "/", "pop_study_",scenario,"_"
                               ,as.character(N),"_",as.character(nsim),".RData"))
  }else if (design_type =="GP")
{
  final_data <- gp_study()
}else
  stop("Insert a valid design type")


###########################################################
######   SCENARIO B #######
###########################################################
# 10 000
N <- 100

# 6) how many individuals?
nsim <- 10000
scenario <- "B"
load(paste0("Simulation/data for simulation/scenario_",scenario,".RData"))

sim_test_listB <- sim_mm_traj(nsim,N,scenario_obj = scenario_obj_B,seed=1258,scenario = "B")
sim_testB <- sim_test_listB[[1]]
sim_MPB <- sim_test_listB[[2]]

save(sim_testB, file = paste0(output_folder, "/", "sim_data_",scenario,"_",as.character(N),"_",as.character(nsim),".RData"))
save(sim_MPB, file = paste0(output_folder, "/", "sim_MP_",scenario,"_",as.character(N),"_",as.character(nsim),".RData"))
message("Dataset saved in 'Simulation Outputs' folder.")

pop <- popbased_study(sim_testB,sim_MPB, underreporting = FALSE,
                      obs_t_int = obs_t_int, age_tr = age_tr)

save(pop, file = paste0(output_folder, "/", "pop_study_",scenario,"_"
                        ,as.character(N),"_",as.character(nsim),".RData"))

rv_pop <- rv_study_new(sim_testB,sim_MPB, underreporting = FALSE)

save(rv_pop, file = paste0(output_folder, "/", "rv_pop_study_",scenario,"_"
                        ,as.character(N),"_",as.character(nsim),".RData"))

# scenario da 3000
# 1) How many datasets?
N <- 100

# 6) how many individuals?
nsim <- 3000

seed = 1259

sim_test_listB <- sim_mm_traj(nsim,N,scenario_obj = scenario_obj_B,seed=seed,scenario = "B")
sim_testB <- sim_test_listB[[1]]
sim_MPB <- sim_test_listB[[2]]

save(sim_testB, file = paste0(output_folder, "/", "sim_data_",scenario,"_",as.character(N),"_",as.character(nsim),".RData"))
save(sim_MPB, file = paste0(output_folder, "/", "sim_MP_",scenario,"_",as.character(N),"_",as.character(nsim),".RData"))
message("Dataset saved in 'Simulation Outputs' folder.")

pop <- popbased_study(sim_testB,sim_MPB, underreporting = FALSE,
                      obs_t_int = obs_t_int, age_tr = age_tr)

save(pop, file = paste0(output_folder, "/", "pop_study_",scenario,"_"
                        ,as.character(N),"_",as.character(nsim),".RData"))
# random visit scenario

rv_pop <- rv_study_new(sim_testB,sim_MPB, underreporting = FALSE)

save(rv_pop, file = paste0(output_folder, "/", "rv_pop_study_",scenario,"_"
                           ,as.character(N),"_",as.character(nsim),".RData"))


#scenario 100 000
N <- 100
nsim <- 100000
seed <- 1260
path_sim <- "Simulation"
output_folder <-paste0(path_sim, "/", "Simulation Outputs")

sim_test_listB <- sim_mm_traj(nsim,N,scenario_obj = scenario_obj_B,seed=seed,scenario = "B")
sim_testB <- sim_test_listB[[1]]
sim_MPB <- sim_test_listB[[2]]

save(sim_testB, file = paste0(output_folder, "/", "sim_data_",scenario,"_",as.character(N),"_",as.character(nsim),".RData"))
save(sim_MPB, file = paste0(output_folder, "/", "sim_MP_",scenario,"_",as.character(N),"_",as.character(nsim),".RData"))
message("Dataset saved in 'Simulation Outputs' folder.")

pop <- popbased_study(sim_testB,sim_MPB, underreporting = FALSE,
                      obs_t_int = obs_t_int, age_tr = age_tr)

save(pop, file = paste0(output_folder, "/", "pop_study_",scenario,"_"
                        ,as.character(N),"_",as.character(nsim),".RData"))