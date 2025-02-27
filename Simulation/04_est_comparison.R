library(ggplot2)
library(dplyr)
library(survminer)
library(tidyverse) 

########## this script is meant for visualizing estimantes/estimands only ########

# load the true param
pr <-load(paste0("results/results_estimates/true_param.RData"))
true_param <- get(pr)
source("Functions//plot_estimates.R")
source("Functions//performance_measures.R")

result_folder <-  "results/results_performance"
save_est <- FALSE
if (!dir.exists(result_folder)) {
  dir.create(result_folder)
  message("Folder 'results_performance' created.")
} else {
  message("Folder 'results_performance' already exists.")
}


###############################################
################## N = 3000 ################### 
###############################################
scenario <- "B"
nsim <- 3000
# load population study estimates
st <- "ps"
pt <-load(paste0("results/results_estimates/model_est_",st,"_",scenario,"_",nsim,".RData"))
param <- get(pt)
param_ps <- param$param_df
conv_time <- param$time_df

param_ps$model <- dplyr::recode(param_ps$model,
                                "flexsurv"="PS_ETESMM",
                                "flexsurv_true" = "PS_ETESMM_base",
                                "flexsurv_imp" = "PS_ImpETESMM",
                                "msm_base"= "PS_ApproxTIMM",
                                "msm_misc"= "PS_ApproxTIHMM",
                                "nhm_base"= "PS_TIMM",
                                "nhm_misc"= "PS_TIHMM"
                                
)

param_rec_states <- param_ps%>% filter(model!="PS_ETESMM_base")
boxplot_param(param_rec_states,true_param,nsim)
dotplot_param(param_rec_states,true_param,nsim)

pm_df <- compute_pm(param_rec_states,true_param) 
dotplot_rel_bias(param_rec_states, true_param,nsim)
param_rec_states_filt <- param_rec_states %>% filter(trans=="1->2")
dotplot_rel_bias_misc(param_rec_states_filt, true_param,nsim)

# load random visits estimates 
st <- "rv"
pt <-load(paste0("results/results_estimates/model_est_",st,"_",scenario,"_",nsim,".RData"))
param <- get(pt)
param_rv <- param$param_df
param_rv$model <- dplyr::recode(param_rv$model,
                                "flexsurv"="IV_ETESMM",
                                "flexsurv_true" = "IV_ETESMM_base",
                                "flexsurv_imp" = "IV_ImpETESMM",
                                "msm_base"= "IV_ApproxTIMM",
                                "msm_misc"= "IV_ApproxTIHMM",
                                "nhm_base"= "IV_TIMM",
                                "nhm_misc"= "IV_TIHMM"
                                
)
dotplot_param(param_rv,true_param,nsim)
dotplot_rel_bias(param_rv, true_param,nsim)

# load exact times estimates 
st <- "ET"
pt2 <-load(paste0("results/results_estimates/model_est_",st,"_",scenario,"_",nsim,".RData"))
param <- get(pt2)
param_df_ET <- param$param_df
conv_time_ET <- param$time_df

param_df_ET$model <- dplyr::recode(param_df_ET$model,
                                            "flexsurv"="ETRC",
                                            "flexsurv_true" = "ETES"
)
dotplot_param(param_df_ET,true_param,nsim)


param_df_3000 <- rbind(param_ps, param_rv, param_df_ET)
param_df_3000<- param_df_3000%>% filter(!(model %in% c("ETRC", "IV_ETESMM_base", "PS_ETESMM_base")))
dotplot_param(param_df_3000,true_param,3000)
dotplot_rel_bias_misc(param_df_3000, true_param,3000,"PS")
dotplot_rel_bias_misc(param_df_3000, true_param,3000,"IV")
dotplot_rel_bias_1_2(param_df_3000, true_param,3000)
dotplot_rel_bias_1_2(param_df_3000, true_param,3000,"PS")

pm_df_3000 <- compute_pm(param_df_3000,true_param) 
if (save_est){
  save(pm_df_3000, file =paste0(result_folder,"/pm_df_3000.RData"))
}

# stability of estimates as n of datasets increases:
m_names <- unique(param_df_3000$model)
for (i in 1:length(m_names)){
  model_name <-m_names[i]
  df <- param_df_3000%>% filter(model==model_name)
  print(plot_est_stability(df, true_param, model_name,nsim=3000))
}
# coverage stability
for (i in 1:length(m_names)){
  model_name <-m_names[i]
  df <- param_df_3000%>% filter(model==model_name)
  print(plot_cov_stability(df, true_param, model_name,nsim=3000))
}
df <- param_df_3000%>% filter(model=="ETES")
plot_est_stability2(df, true_param, "ETES",nsim=3000)

###############################################
################## N = 10000 ##################
###############################################
scenario <- "B"
nsim <- 10000
# load population study estimates
st <- "ps"
pt <-load(paste0("results/results_estimates/model_est_",st,"_",scenario,"_",nsim,".RData"))
param <- get(pt)
param_ps <- param$param_df
conv_time <- param$time_df

# remove flexsurv_true and compare the models fitted after LCA
param_ps$model <- dplyr::recode(param_ps$model,
                                "flexsurv"="PS_ETESMM",
                                "flexsurv_true" = "PS_ETESMM_base",
                                "flexsurv_imp" = "PS_ImpETESMM",
                                "msm_base"= "PS_ApproxTIMM",
                                "msm_misc"= "PS_ApproxTIHMM",
                                "nhm_base"= "PS_TIMM",
                                "nhm_misc"= "PS_TIHMM",
                                
)
param_rec_states <- param_ps%>% filter(model!="PS_ETESMM_base")

boxplot_param(param_rec_states,true_param,nsim)
dotplot_param(param_rec_states,true_param,nsim)

pm_df <- compute_pm(param_rec_states,true_param) 
dotplot_rel_bias(param_rec_states, true_param,nsim)
param_rec_states_filt <- param_rec_states %>% filter(trans=="1->2")
dotplot_rel_bias_misc(param_rec_states_filt, true_param,nsim)

# load random visits estimates 
st <- "rv"
pt <-load(paste0("results/results_estimates/model_est_",st,"_",scenario,"_",nsim,".RData"))
param <- get(pt)
param_rv <- param$param_df
dotplot_param(param_rv,true_param,nsim)
dotplot_rel_bias(param_rv, true_param,nsim)

param_rv$model <- dplyr::recode(param_rv$model,
                                "flexsurv"="IV_ETESMM",
                                "flexsurv_true" = "IV_ETESMM_base",
                                "flexsurv_imp" = "IV_ImpETESMM",
                                "msm_base"= "IV_ApproxTIMM",
                                "msm_misc"= "IV_ApproxTIHMM",
                                "nhm_base"= "IV_TIMM",
                                "nhm_misc"= "IV_TIHMM",
                                
)
# load exact times estimates 
st <- "ET"
pt2 <-load(paste0("results/results_estimates/model_est_",st,"_",scenario,"_",nsim,".RData"))
param <- get(pt2)
param_df_ET <- param$param_df
conv_time_ET <- param$time_df

param_df_ET$model <- dplyr::recode(param_df_ET$model,
                                   "flexsurv"="ETRC",
                                   "flexsurv_true" = "ETES"
)
dotplot_param(param_df_ET,true_param,nsim)


param_df_10000 <- rbind(param_ps, param_rv, param_df_ET)
param_df_10000<- param_df_10000%>% filter(!(model %in% c("ETRC", "IV_ETESMM_base", "PS_ETESMM_base")))
dotplot_param(param_df_10000,true_param,nsim)
dotplot_rel_bias_misc(param_df_10000, true_param,10000, "PS")
dotplot_rel_bias_misc(param_df_10000, true_param,10000, "IV")
dotplot_rel_bias_1_2(param_df_10000, true_param,10000)
pm_df_10000 <- compute_pm(param_df_10000,true_param) 
if (save_est){
  save(pm_df_10000, file =paste0(result_folder,"/pm_df_10000.RData"))
}

pm_df_10000_misc <- pm_df_10000 %>% filter(trans=="1->2", parameter %in% c("beta_educ", "beta_sex"),  model %in% c("IV_TIHMM","IV_TIMM", "PS_TIHMM","PS_TIMM","IV_ApproxTIHMM","IV_ApproxTIMM", "PS_ApproxTIHMM","PS_ApproxTIMM"))

#pm_df_10000_misc <- pm_df_10000 %>% filter(trans=="1->2", model %in% c("IV_ApproxTIHMM","IV_ApproxTIMM", "PS_ApproxTIHMM","PS_ApproxTIMM"))


# stability of estimates as n of datasets increases:
m_names <- unique(param_df_10000$model)
for (i in 1:length(m_names)){
  model_name <-m_names[i]
  df <- param_df_10000%>% filter(model==model_name)
  print(plot_est_stability(df, true_param, model_name,nsim=10000))
}
df <- param_df_10000%>% filter(model=="ETES")
plot_est_stability2(df, true_param, "ETES",nsim=10000)
# coverage stability
for (i in 1:length(m_names)){
  model_name <-m_names[i]
  df <- param_df_10000%>% filter(model==model_name)
  print(plot_cov_stability(df, true_param, model_name,nsim=10000))
}
################## flexsurv comparison #################
# compare pop_study, exact times and random visits fitted both on sim_MP (exact states) and MP (recovered states)

param_ps <- param_df %>% filter(model %in% c ("flexsurv","flexsurv_true", "flexsurv_imp"))
param_rv <- param_df %>% filter(model %in% c ("flexsurv","flexsurv_true", "flexsurv_imp"))


param_df_flex <- rbind(param_df_ET, param_ps,param_rv)

dotplot_param(param_df_flex,true_param,nsim)

# compute performance measures
pm_df <- compute_pm(param_df_flex,true_param) 
dotplot_rel_bias(param_df_flex, true_param,nsim)



################## Benchmark model #################
m1 <- param_df_3000 %>% filter(model=="ETES")%>%
  mutate(nsim="3000")
m2 <-  param_df_10000 %>% filter(model=="ETES")%>%
  mutate(nsim="10000")
param_bench <- rbind(m1,m2)
param_bench <- param_bench %>%
  mutate(nsim=factor(nsim,levels = c("3000", "10000")))
boxplot_param_ss(param_bench, true_param)
dotplot_param_ss(param_bench, true_param)
dotplot_rel_bias_ss(param_bench, true_param)

# table of performance measure for benchmark model:

pm_df_bench <- rbind(compute_pm(m1,true_param),compute_pm(m2,true_param))
pm_df_bench
pm_df_bench <- pm_df_bench %>% dplyr::select(-model)
pm_df_bench[, sapply(pm_df_bench, is.numeric)] <- 
  lapply(pm_df_bench[, sapply(pm_df_bench, is.numeric)], function(x) format(round(x, 3), nsmall = 3))
latex_table <- xtable(pm_df_bench)
digits_vector <- c(0, ifelse(sapply(pm_df_bench, is.numeric), 3, 0))
print(latex_table, type = "latex", include.rownames = FALSE,digits = digits_vector)

############################## sample size comparison ################
pr <-load(paste0("results/results_estimates/true_param.RData"))
true_param <- get(pr)

scenario <- "B"
nsim <- 10000
# load population study estimates
st <- "ps"
pt <-load(paste0("results/results_estimates/model_est_",st,"_",scenario,"_",nsim,".RData"))
param <- get(pt)
param_df_2 <- param$param_df
conv_time_2 <- param$time_df
param_rc_2 <- param_df_2%>% filter(model!="flexsurv_true")
dotplot_param(param_rc_2,true_param,nsim)
pm_df <- compute_pm(param_rc_2,true_param) 
dotplot_rel_bias(param_rc_2, true_param,nsim)
param_rc_fit_2 <- param_rc_2 %>% filter(trans=="1->2")
dotplot_rel_bias_misc(param_rc_fit_2, true_param,nsim)



#####################################################
################# Hazard plots ##################### 
source("Functions//plot_hazard.R")
t_vals <- seq(60, 90, length.out = 100)  # Range of x values

#plot_cum_hazard(pm_df_3000,true_param,t_vals, 3000)
#plot_cum_hazard_gg(pm_df_3000,true_param,t_vals, nsim= 3000)
# single plot
#plot_hazard_gg(pm_df_3000,true_param,t_vals, nsim= 3000, "PS")
plot_hazard_gg_2(pm_df_3000,true_param,t_vals, nsim= 3000, "PS")
plot_hazard_gg_2(pm_df_3000,true_param,t_vals, nsim= 3000, "IV")

#plot_cum_hazard_gg(pm_df_10000,true_param,t_vals, nsim= 10000)
plot_hazard_gg_2(pm_df_10000,true_param,t_vals, nsim= 10000, "PS")
plot_hazard_gg_2(pm_df_10000,true_param,t_vals, nsim= 10000, "IV")

t_vals2 <- seq(60, 100, length.out = 100)
base_meas_3000 <- eval_cum_hazard(pm_df_3000,true_param,t_vals, nsim= 3000)
base_meas_10000 <- eval_cum_hazard(pm_df_10000,true_param,t_vals, nsim= 10000)

eval_cum_hazard(pm_df_3000,true_param,t_vals2, nsim= 3000, "PS")
base_meas_3000
base_meas_3000[, sapply(base_meas_3000, is.numeric)] <- 
  lapply(base_meas_3000[, sapply(base_meas_3000, is.numeric)], function(x) {
    ifelse(abs(x) < 0.001, formatC(x, format = "e", digits = 3), round(x, 4))
  })
base_meas_10000[, sapply(base_meas_10000, is.numeric)] <- 
  lapply(base_meas_10000[, sapply(base_meas_10000, is.numeric)], function(x) {
    ifelse(abs(x) < 0.001, formatC(x, format = "e", digits = 3), round(x, 4))
  })


library(xtable)
latex_table <- xtable(base_meas_3000)
print(latex_table, type = "latex", include.rownames = FALSE)
latex_table <- xtable(base_meas_10000)
print(latex_table, type = "latex", include.rownames = FALSE)

################################################################
###################### Latex tables ###########################
library(xtable)
# only ps, covariates, trans 1->2
pm_df_10000_misc <- pm_df_10000 %>% filter(trans=="1->2", parameter %in% c("beta_educ", "beta_sex"),  model %in% c("PS_TIHMM","PS_TIMM","PS_ApproxTIHMM","PS_ApproxTIMM","PS_ETESMM" ,"PS_ImpETESMM")) %>%
  arrange(parameter) %>% dplyr::select(-trans)
pm_df_10000_misc[, sapply(pm_df_10000_misc, is.numeric)] <- 
  lapply(pm_df_10000_misc[, sapply(pm_df_10000_misc, is.numeric)], function(x) format(round(x, 3), nsmall = 3))

latex_table <- xtable(pm_df_10000_misc)
digits_vector <- c(0, ifelse(sapply(pm_df_10000_misc, is.numeric), 3, 0))
print(latex_table, type = "latex", include.rownames = FALSE,digits = digits_vector)

# only ps, covariates, trans 1->2

pm_df_10000_misc <- pm_df_10000 %>% filter(trans=="1->2", parameter %in% c("beta_educ", "beta_sex"),  model %in% c("IV_TIHMM","IV_TIMM","IV_ApproxTIHMM","IV_ApproxTIMM", "IV_ImpETESMM", "IV_ETESMM")) %>%
  arrange(parameter) %>% dplyr::select(-trans)
pm_df_10000_misc[, sapply(pm_df_10000_misc, is.numeric)] <- 
  lapply(pm_df_10000_misc[, sapply(pm_df_10000_misc, is.numeric)], function(x) format(round(x, 3), nsmall = 3))

latex_table <- xtable(pm_df_10000_misc)
digits_vector <- c(0, ifelse(sapply(pm_df_10000_misc, is.numeric), 3, 0))
print(latex_table, type = "latex", include.rownames = FALSE,digits = digits_vector)

# rate and shape:
pm_df_10000_base <- pm_df_10000 %>% filter( parameter %in% c("rate", "shape"),  model %in% c("PS_TIHMM","PS_TIMM","PS_ApproxTIHMM","PS_ApproxTIMM","PS_ETESMM" ,"PS_ImpETESMM")) %>%
  arrange(parameter)
pm_df_10000_base[, sapply(pm_df_10000_base, is.numeric)] <- 
  lapply(pm_df_10000_base[, sapply(pm_df_10000_base, is.numeric)], function(x) format(round(x, 3), nsmall = 3))

latex_table <- xtable(pm_df_10000_base)
digits_vector <- c(0, ifelse(sapply(pm_df_10000_base, is.numeric), 3, 0))
print(latex_table, type = "latex", include.rownames = FALSE,digits = digits_vector)

### Nsim = 3000
pm_df_3000_misc <- pm_df_3000 %>% filter(trans=="1->2", parameter %in% c("rate", "shape"),  model %in% c("PS_TIHMM","PS_TIMM","PS_ApproxTIHMM","PS_ApproxTIMM","PS_ETESMM" ,"PS_ImpETESMM")) %>%
  arrange(parameter) %>% dplyr::select(-trans)
pm_df_3000_misc[, sapply(pm_df_3000_misc, is.numeric)] <- 
  lapply(pm_df_3000_misc[, sapply(pm_df_3000_misc, is.numeric)], function(x) format(round(x, 3), nsmall = 3))

latex_table <- xtable(pm_df_3000_misc)
digits_vector <- c(0, ifelse(sapply(pm_df_3000_misc, is.numeric), 3, 0))
print(latex_table, type = "latex", include.rownames = FALSE,digits = digits_vector)


pm_df_3000_misc <- pm_df_3000 %>% filter(trans=="1->2", parameter %in% c("rate", "shape"),  model %in% c("IV_TIHMM","IV_TIMM","IV_ApproxTIHMM","IV_ApproxTIMM", "IV_ImpETESMM", "IV_ETESMM")) %>%
  arrange(parameter) %>% dplyr::select(-trans)
pm_df_3000_misc[, sapply(pm_df_3000_misc, is.numeric)] <- 
  lapply(pm_df_3000_misc[, sapply(pm_df_3000_misc, is.numeric)], function(x) format(round(x, 3), nsmall = 3))

latex_table <- xtable(pm_df_3000_misc)
digits_vector <- c(0, ifelse(sapply(pm_df_3000_misc, is.numeric), 3, 0))
print(latex_table, type = "latex", include.rownames = FALSE,digits = digits_vector)
