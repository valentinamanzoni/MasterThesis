
create_baseline_table <- function(pop, index=1){
  summ_sim <- pop %>%
    group_by(patient_id) %>%
    slice(index) %>%  # Select first observation per patient
    ungroup() %>% 
    dplyr::select(-dataset_id, -patient_id, -Age_entry, -Age_death, -dm_sex,-educ_el, -time_in_study, -ndis, -Age_exit, -age, -visit_number,-age_group
                  , - base_MP, - MP_base, -lag_age, -MP_sim)
  
  summary_table_sim <- tbl_summary(summ_sim)
  
  return(summary_table_sim)
}

create_age_group_table <- function(pop){
  # disease freq by age group
  sum_pop_ag <- pop %>%
    group_by(patient_id, age_group) %>%
    # keep only last obs where I have the same value of age_group multiple times for the same patient: 
    slice_tail(n = 1) %>% 
    ungroup() %>%
    dplyr::select(-dataset_id, -patient_id, -Age_entry, -Age_death, -dm_sex,-educ_el, -time_in_study, -ndis, -Age_exit, -age, -visit_number
                  , - base_MP, - MP_base, -lag_age, -MP_sim)
  
  summary_table_sim <- tbl_summary(sum_pop_ag, by = age_group, missing = 'no')
  
  return(summary_table_sim)
}
