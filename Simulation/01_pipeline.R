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

# save simulation results 
output_folder <-paste0(path_sim, "/", "Simulation Outputs")

if (!dir.exists(output_folder)) {
  dir.create(output_folder)
  message("Folder 'Simulation Outputs' created.")
} else {
  message("Folder 'Simulation Outputs' already exists.")
}

# years in between follow-ups and corresponding age threshold for population based scenario
# patients with age <= 78 years have follow ups every 6 years, > 78 every 3 years
age_tr <- 78
obs_t_int <- c(6,3)

scenario <- "B"
load(paste0("Simulation/data for simulation/scenario_",scenario,".RData"))


###########################################################
######   SCENARIO B: 2 Multimorbidity latent classes #######
###########################################################

# 10000 subject scenario
# 1) How many datasets?
N <- 100

# 2) how many individuals?
nsim <- 10000

# Set seed for reproduce results 
seed<- 1258

sim_test_listB <- sim_mm_traj(nsim,N,scenario_obj = scenario_obj_B,seed=seed,scenario = "B")
sim_testB <- sim_test_listB[[1]]
sim_MPB <- sim_test_listB[[2]]

save(sim_testB, file = paste0(output_folder, "/", "sim_data_",scenario,"_",as.character(N),"_",as.character(nsim),".RData"))
save(sim_MPB, file = paste0(output_folder, "/", "sim_MP_",scenario,"_",as.character(N),"_",as.character(nsim),".RData"))
message("Dataset saved in 'Simulation Outputs' folder.")

### Population based study

pop <- popbased_study(sim_testB,sim_MPB, underreporting = FALSE,
                      obs_t_int = obs_t_int, age_tr = age_tr)

save(pop, file = paste0(output_folder, "/", "pop_study_",scenario,"_"
                        ,as.character(N),"_",as.character(nsim),".RData"))

### Irregular visits

rv_pop <- rv_study_new(sim_testB,sim_MPB, underreporting = FALSE)

save(rv_pop, file = paste0(output_folder, "/", "rv_pop_study_",scenario,"_"
                        ,as.character(N),"_",as.character(nsim),".RData"))

# 3000 subject scenario
# 1) How many datasets?
N <- 100

# 2) how many individuals?
nsim <- 3000

seed <- 1259

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


#######################################
#######  UNDERREPORTING ###############
#######################################

source("Simulation/fun/Study_design_fun_2.R")
set.seed(42)

# 1) which study type?
study_type <- "pop_study"

# 2) How many datasets?
N <- 100

# 3) how many individuals?
nsim <- 3000

# 4) Which scenario?
scenario <- "B"

# load the data 
load(paste0("Simulation/Simulation Outputs/",study_type,"_",scenario,"_"
            ,as.character(N),"_",as.character(nsim),".RData"))

# 5) list of underreported disease
underrep_diseases <- c("dementia", "depression_mood_dis", "chronic_kidney_dis", "osteoarthr_degen_joint_dis", "deafness_hearing_loss")
# 6) vector with corresponding underreporting probabilities
diseases_prob <- c(0.2,0.2,0.2,0.2, 0.2)

pop_under <- apply_underrep(pop, underrep_diseases , diseases_prob )

save(pop_under, file = paste0(output_folder, "/underrep_", study_type,"_",scenario,"_"
                           ,as.character(N),"_",as.character(nsim),".RData"))

