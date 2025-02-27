
library(data.table)
library(dplyr)
library(future.apply)



# load simulated Data
# base_path <- "~/Dropbox/Il mio Mac (AirdiValentina.homenet.telecomitalia.it)/Desktop/Tesi/Codice librerie/Code/Tesi Valentina/Data"
# full_path <- file.path(base_path, "scenario_A.RData")
# load(full_path)
# 


###############   SCENARY 1 #############
############## POPULATION-BASED STUDY ##########
# popbased study returns a dataframe of simulated data of a population-base study where
# data is the underlying simulated population
# underreporting = TRUE if we want to include underepporting for some of the simulated diseases
# underrep_dis vector containing the names of underreported diseases
# underrep_pop vector with the probabilities of underepporting for each disease
# specified in underrep_dis
# obs_t_int vector containing information about the period (in years) between follow-ups, can contain one or 2 values
# age_tr if more than one value is passed to obs_t_int then it needs to contain one value indicating after which age the follow up period changes


popbased_study <- function (data, sim_MP, underreporting = FALSE, underrep_dis = NULL, underrep_prob = NULL, obs_t_int, age_tr= NA){
  # filter data to keep observations according to obs_t_init
  data <- as.data.frame(data)
  # DISCRETIZATION : one row per follow up 
  
  if (length(obs_t_int)>1 & is.na(age_tr))
  {
    stop("ERROR. obs_t_int and age_tr have incosistent numebr of elements.")
  }
  else if(length(obs_t_int)==1)
  { data_long <- data %>%
      group_by(dataset_id, patient_id) %>%
      # MODIFY ONCE FIRST PART OF THE SIM IS FIXED
      mutate(age = list(if (Age_entry <= Age_death) seq(Age_entry, Age_exit, by = obs_t_int) else NA)) %>%
      # mutate(age = list(seq(Age_entry, Age_death, by = obs_t_int))) %>%
      unnest(age) %>%
      mutate(visit_number = row_number()) %>%
      ungroup() 
  }
  
  else {
    create_visit_seq <- function (Age_entry,Age_death,obs_t_int, age_tr){
      if (Age_entry <= Age_death)
      { 
        if (Age_entry<age_tr)
          s1 <- seq(Age_entry, min(age_tr+(obs_t_int[1]*0.99), Age_death), by = obs_t_int[1])
        else 
          s1 <- NULL
        last_element <- max(Age_entry,tail(s1, 1))
        if (Age_death > last_element)
          s2 <- seq(last_element, Age_death, by = obs_t_int[2])
        else 
          s2 <-NULL
        sol <- unique(c(s1, s2))
        return(sol)
      } else {return(NA)}
    }
    
    data_long <- data %>%
      rowwise() %>%
      mutate(age = list(create_visit_seq(Age_entry,Age_exit,obs_t_int, age_tr))) %>%
      unnest(age) %>%
      ungroup() %>%
      group_by(dataset_id,patient_id) %>%
      mutate(visit_number = row_number()) %>%
      ungroup()
  }

  data_long <- create_age_group(data_long)
  # save the real simulated MP 
  sim_MP_util <- sim_MP %>% dplyr::select(dataset_id, patient_id,Age, MP,lag_age)
  data_long<- add_MP_sim(data_long,sim_MP_util)
  # compare each Age_disease (no Age_entry and Age_death) with age if Age_disease <= age then 1 else 0
  disease_columns <- grep("^Age_", names(data), value = TRUE)
  disease_columns <- setdiff(disease_columns, c("Age_entry", "Age_death", "Age_exit"))

  data_long <- data_long %>%
    mutate(across(all_of(disease_columns), 
                  ~ ifelse(. <= age, 1, 0), 
                  .names = "{sub('Age_', '', .col)}"))
  data_long <- data_long %>% dplyr::select(-all_of(disease_columns))
  
  #diseases_list <-sub("^Age_", "", disease_columns)
  #data_long <- create_age_group(data_long, diseases_list)
  
  if (underreporting){
    if (length(underrep_dis)!=length(underrep_prob)){
      stop("ERROR. underrep_dis and underrep_prob must have same number of element. 
            Specify one value of underrep_prob per disease passed to underrep_dis.")
    }
    # for each disease in underrep_dis simulate from a bernuolli (underrep_prob[i])
    # if outcome = 1 -> disease underreported, then removed from dataset
    for (i in seq_along(underrep_dis)) {
      disease <- underrep_dis[i]
      prob <- underrep_prob[i]
      
      patients_with_disease <- unique(data_long$patient_id[data_long[[disease]] == 1])
      
      # Generate underreporting flags only for patients with the disease
      if (length(patients_with_disease) > 0) {
        underreport_flags <- rbinom(length(patients_with_disease), 1, prob)
        
        for (j in seq_along(patients_with_disease)) {
          pat_id <- patients_with_disease[j]
          indices <- which(data_long$patient_id == pat_id) 
          
          if (length(indices) > 0) {
            data_long$ndis[indices] <- data_long$ndis[indices] - 1 # decrease ndis whenever there's an underreported disease
            # Update the disease column while keeping NAs
            data_long[[disease]][indices] <- ifelse(
              is.na(data_long[[disease]][indices]), 
              NA,  # Keep NA
              ifelse(data_long[[disease]][indices] == 1 & underreport_flags[j] == 1, 0, 1)
            )
          }
        }
      }
    }
  }  
  
  data_long[is.na(data_long)] <- 0
  return(data_long)
}

create_age_group <- function(data){
  dat_ag <- data %>%
    mutate(age_group = case_when(
      age >= 60 & age < 65 ~ "60-65",
      age >= 65 & age < 70 ~ "65-70",
      age >= 70 & age < 75 ~ "70-75",
      age >= 75 & age < 80 ~ "75-80",
      age >= 80 & age < 85 ~ "80-85",
      age >= 85 & age < 90 ~ "85-90",
      age >= 90 ~ "90+",
      TRUE ~ NA_character_  # Handle cases where age is less than 60 or missing
    ))
  return(dat_ag)
}

add_MP_sim <- function(data_long, sim_MP) {
  sim_MP <- sim_MP %>%
    arrange(dataset_id, patient_id, Age) %>%
    mutate(next_Age = lead(Age, default = Inf))
  
  # Join data_long with sim_MP
  data_long <- data_long %>%
    left_join(sim_MP, by = c("dataset_id", "patient_id"), relationship = "many-to-many") %>%
    filter(age >= Age & age < next_Age) %>%
    mutate(MP_sim=MP) %>%
    dplyr::select(-next_Age, -Age, -MP)
  
  return(data_long)
}


###############   SCENARY 2 #############
############ RANDOM VISIT ########
# simulate patient visit to based on a Weibull
# rate dependent on age, sex, disease burden (standard case just uses fixed shape and rate)
# if visit_rate not provided then it is calculated using function calculate_visit_rate ()
calculate_visit_rate <- function(base_rate=NULL, age, sex, socio_ec, num_diseases) {
  
  if (is.null(base_rate)){
    base_rate <- 2   # if not provided assume base rate to be 2 visits per year 
  }

  age_factor <- 0.01 * (age - 60)   # Each year adds 1% to the base rate
  
  # adjust rate based on sex and num of diseases
  sex_factor <- ifelse(sex == "female", 1.2, 1.0)  # Females 20% more likely
  disease_factor <- 1 + 0.3 * num_diseases  # Each disease adds 30% to the base rate
  
  # Total rate
  visit_rate <- base_rate * (1 + age_factor) * sex_factor * disease_factor
  return(visit_rate)
}

rv_study_new <- function(data, sim_MP, underreporting = FALSE, underrep_dis = NULL, underrep_prob = NULL, base_rate=NULL){
  
  data <- as.data.frame(data)
  # simulate from a weibull
  # visits for each patient in between Age_entry and Age_exit
  simulate_patient_visits_old <- function(Age_entry, Age_exit, sex, socio_ec, num_diseases) {
    
    current_age <- Age_entry
    visit_ages <- c()  # Store ages at which visits occur
    
    while (TRUE) {
      
      # Calculate Weibull parameters (scale and shape) using the patient's current state
      #visit_params <- calculate_visit_rate(base_rate = base_rate, age = current_age, sex = sex,socio_ec = socio_ec,num_diseases = num_diseases)
      
      #shape <- visit_params$shape
      #scale <- visit_params$scale
      
      shape<- 5
      scale<- 0.4
      
      time_to_next_visit <- rweibull(1, shape = shape, scale = scale)
      next_visit_age <- current_age + time_to_next_visit
      
      # Stop if the next visit age is beyond Age_exit
      if (next_visit_age > Age_exit) break
      
      visit_ages <- c(visit_ages, next_visit_age)
      current_age <- next_visit_age
    }
    
    return(visit_ages)
  }
  simulate_patient_visits <- function(Age_entry, Age_exit, sex, socio_ec, num_diseases) {
    current_age <- Age_entry
    max_visits <- 30 # max number of visits
    visit_ages <- numeric(max_visits)
    visit_count <- 0
    
    while (TRUE) {
      time_to_next_visit <- rweibull(1, shape = 3, scale = 1.5)
      next_visit_age <- current_age + time_to_next_visit
      if (next_visit_age > Age_exit) break
      
      visit_count <- visit_count + 1
      visit_ages[visit_count] <- next_visit_age
      current_age <- next_visit_age
      
      if (visit_count >= max_visits) break  # Prevent infinite growth
    }
    
    return(visit_ages[1:visit_count])  # Trim unused elements
  }
  
  
  data_long <- data %>%
    group_by(dataset_id, patient_id) %>%
    mutate(age =list(simulate_patient_visits(Age_entry, Age_exit, dm_sex, educ_el, ndis))) %>%
    unnest(age) %>%
    mutate(visit_number = row_number()) %>%
    ungroup() 
  
  data_long <- create_age_group(data_long)
  # save the real simulated MP 
  sim_MP_util <- sim_MP %>% dplyr::select(dataset_id, patient_id,Age, MP,lag_age)
  data_long<- add_MP_sim(data_long,sim_MP_util)
  # compare each Age_disease (no Age_entry and Age_death) with age if Age_disease <= age then 1 else 0
  disease_columns <- grep("^Age_", names(data), value = TRUE)
  disease_columns <- setdiff(disease_columns, c("Age_entry", "Age_death", "Age_exit"))
  
  data_long <- data_long %>%
    mutate(across(all_of(disease_columns), 
                  ~ ifelse(. <= age, 1, 0), 
                  .names = "{sub('Age_', '', .col)}"))
  data_long <- data_long %>% dplyr::select(-all_of(disease_columns))
  
  data_long[is.na(data_long)] <- 0
  return(data_long)
}

################## UNDERREPORTING ###############

apply_underrep <- function(data, disease_names, dis_prob){
  

  for (i in (1:length(disease_names))){
    # estrarre righe del soggetto in cui ha quella malattia, per ogni riga sampling da binom per underrep
    # filtrare da prima riga a cui non c'è underrep in poi e mettere flag a 0 
    # 
    dis <- disease_names[i]
    if (!(dis %in% colnames(data))){
      message(paste("No column", dis, "in the data. "))
      next
    }
    prob <- dis_prob[i]
    
    data_2 <- data %>% filter(.data[[dis]]==1)
    data_2$flag <- rbinom(nrow(data_2),1,prob)
    
    data_2 <- data_2 %>% dplyr::group_by(dataset_id, patient_id) %>%
      mutate(flag2 = case_when(lag(flag, default=2)==1 & flag==0 ~ 1, lag(flag, default=2)==0 ~ 2,max(flag)==0 ~ 2, TRUE ~ 0)) 
    
    data_2 <- data_2%>%
      filter(flag2>0) %>%
      summarise(min_no_under= min(age))
    
    data <- data%>% left_join(data_2)
    
    data[dis] <- case_when(is.na(data$min_no_under)  ~ 0, data$age >= data$min_no_under ~ 1, TRUE ~ 0)
    
    data <- data %>% dplyr::select(-min_no_under)
  
  }
  return(data)
}


