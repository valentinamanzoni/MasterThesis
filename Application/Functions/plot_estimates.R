

boxplot_param <- function(param_df,true_param,nsim){
  data_long <- param_df %>%
    pivot_longer(cols = c(rate, shape, beta_educ,beta_sex), 
                 names_to = "parameter", 
                 values_to = "value")
  
  true_param_long <- true_param %>%
    pivot_longer(cols = c(rate, shape, beta_educ, beta_sex), 
                 names_to = "parameter", 
                 values_to = "true_value")
  
  data_long <- data_long %>%
    left_join(true_param_long, by = c("trans", "parameter"))
  
  ggplot(data_long, aes(x = model, y = value, fill = model)) +
    geom_boxplot() +
    geom_hline(aes(yintercept = true_value), color = "red", linetype = "dashed", linewidth = 1) +
    facet_grid(parameter ~ trans, scales = "free_y") +
    theme_bw() +
    labs(
      title = paste("Boxplots of Parameters Grouped by Model and Transition. N=",nsim),
      x = "Model name",
      y = "Estimated parameter"
    )+  
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set2")
  
}

dotplot_param <- function(param_df,true_param,nsim){
  data_long <- param_df %>%
    pivot_longer(cols = c(rate, shape, beta_educ,beta_sex), 
                 names_to = "parameter", 
                 values_to = "value") 
  data_long_se <- param_df %>%
    pivot_longer(cols = c(rate_se, shape_se, beta_educ_se,beta_sex_se), 
                 names_to = "parameter", 
                 values_to = "value_se")%>%
    mutate(parameter = sub("_se$", "", parameter))
  
  data_long$value_se <- data_long_se$value_se
  
  true_param_long <- true_param %>%
    pivot_longer(cols = c(rate, shape, beta_educ, beta_sex), 
                 names_to = "parameter", 
                 values_to = "true_value")
  
   
  data_summary <- data_long %>%
    group_by(model, trans, parameter) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      within_se = mean(value_se, na.rm = TRUE),
      between_se = sd(value, na.rm = TRUE),
      se_value = sqrt(within_se^2 + between_se^2 / n()),
      .groups = "drop"
    )
  data_summary <- data_summary %>%
    left_join(true_param_long, by = c("trans", "parameter"))
  
  # dot plot of mean values with error bars
  ggplot(data_summary, aes(x = model, y = mean_value)) +
    geom_point(size = 3) +  
    geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.2) +  
    geom_hline(aes(yintercept = true_value), color = "red", linetype = "dashed", linewidth = 1) +  
    facet_grid(parameter ~ trans, scales = "free") +
    theme_bw() +
    labs(
      title = paste("Mean estimated parameter grouped by model and transition. N=", nsim),
      x = "Model Name",
      y = "Mean Estimated Parameter"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+
    #scale_color_brewer(palette = "Set2")

}

boxplot_param_ss <- function(param_df,true_param){
  data_long <- param_df %>%
    pivot_longer(cols = c(rate, shape, beta_educ,beta_sex), 
                 names_to = "parameter", 
                 values_to = "value")
  
  true_param_long <- true_param %>%
    pivot_longer(cols = c(rate, shape, beta_educ, beta_sex), 
                 names_to = "parameter", 
                 values_to = "true_value")
  
  data_long <- data_long %>%
    left_join(true_param_long, by = c("trans", "parameter"))
  
  ggplot(data_long, aes(x = nsim, y = value, fill = nsim)) +
    geom_boxplot() +
    geom_hline(aes(yintercept = true_value), color = "red", linetype = "dashed", linewidth = 1) +
    facet_grid(parameter ~ trans, scales = "free_y") +
    theme_bw() +
    labs(
      title = paste("Boxplots of parameters estimates from Benchmark model grouped by sample size and transition"),
      x = "Subjects Sample size",
      y = "Estimated parameter"
    )+  
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position="none") +
    scale_fill_manual(values = c("#29AB87","#006400","#00A86B","#2E8B57","#228B22"))
  
  
}

dotplot_param_ss  <- function(param_df,true_param){
  data_long <- param_df %>%
    pivot_longer(cols = c(rate, shape, beta_educ,beta_sex), 
                 names_to = "parameter", 
                 values_to = "value") 
  data_long_se <- param_df %>%
    pivot_longer(cols = c(rate_se, shape_se, beta_educ_se,beta_sex_se), 
                 names_to = "parameter", 
                 values_to = "value_se")%>%
    mutate(parameter = sub("_se$", "", parameter))
  
  data_long$value_se <- data_long_se$value_se
  
  true_param_long <- true_param %>%
    pivot_longer(cols = c(rate, shape, beta_educ, beta_sex), 
                 names_to = "parameter", 
                 values_to = "true_value")
  
  
  data_summary <- data_long %>%
    group_by(nsim, trans, parameter) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      within_se = mean(value_se, na.rm = TRUE),
      between_se = sd(value, na.rm = TRUE),
      se_value = sqrt(within_se^2 + between_se^2 / n()),
      .groups = "drop"
    )
  data_summary <- data_summary %>%
    left_join(true_param_long, by = c("trans", "parameter"))
  
  # dot plot of mean values with error bars
  ggplot(data_summary, aes(x = nsim, y = mean_value, color=nsim)) +
    geom_point(size = 3) +  
    geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.2) +  
    geom_hline(aes(yintercept = true_value), color = "red", linetype = "dashed", linewidth = 1) +  
    facet_grid(parameter ~ trans, scales = "free") +
    theme_bw() +
    labs(
      title = paste("Mean estimated parameter of Benchmark model grouped by sample size and transition"),
      x = "Subjects Sample size",
      y = "Mean Estimated Parameter"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position="none") +
    scale_color_manual(values = c("#29AB87", "#006400","#00A86B","#2E8B57","#228B22"))
  
}
dotplot_rel_bias_ss <- function(df, true_param){
  data_long <- df %>%
    pivot_longer(cols = c(rate, shape, beta_educ,beta_sex), 
                 names_to = "parameter", 
                 values_to = "value")
  
  true_param_long <- true_param %>%
    pivot_longer(cols = c(rate, shape, beta_educ, beta_sex), 
                 names_to = "parameter", 
                 values_to = "true_value")
  
  data_long <- data_long %>%
    left_join(true_param_long, by = c("trans", "parameter"))
  
  data_summary <- data_long %>%
    group_by(nsim, trans, parameter) %>%
    summarise(
      n = n(), 
      rel_bias = mean((value - true_value)/abs(true_value), na.rm = TRUE),
      se_rel_bias = sd((value - true_value)/abs(true_value), na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  # dot plot of mean values with error bars
  ggplot(data_summary, aes(x = nsim, y = rel_bias, color = nsim)) +
    geom_point(size = 3) +  
    geom_errorbar(aes(ymin = rel_bias - se_rel_bias, ymax =  rel_bias + se_rel_bias), width = 0.2) +  
    geom_hline(aes(yintercept = 0), color = "gray", linetype = "dashed", linewidth = 1) +  
    #facet_wrap(~ trans, scales = "free_y", nrow = 1) +
    facet_wrap(~parameter + trans, scales = "free_y", nrow = 4) +
    theme_bw() +
    labs(
      title = paste("Relative bias of Benchmark model grouped by sample size and transition"),
      x = "Subjects Sample size",
      y = "Relative Bias"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"  # This removes the legend
    ) +
    scale_color_manual(values = c("#29AB87", "#006400","#00A86B","#2E8B57","#228B22"))
  
  
}

plot_param_hist <- function(param_df,true_param,nsim){
data_long <- param_df %>%
  pivot_longer(cols = c(rate, shape, beta_educ, beta_sex), 
               names_to = "parameter", 
               values_to = "value")

true_param_long <- true_param %>%
  pivot_longer(cols = c(rate, shape, beta_educ, beta_sex), 
               names_to = "parameter", 
               values_to = "true_value")
data_long <- data_long %>%
  left_join(true_param_long, by = c("trans", "parameter"))

ggplot(data_long, aes(x = value, fill = model)) +
  geom_histogram(binwidth = 20, alpha = 0.7, position = "identity") +
  geom_vline(aes(xintercept = true_value), color = "red", linetype = "dashed", linewidth = 1) +
  facet_grid(parameter ~ trans, scales = "free") +
  theme_bw() +
  labs(
    title = paste("Histograms of Parameters Grouped by Model and Transition. N =", nsim),
    x = "Parameter Value",
    y = "Count"
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )
}

plot_param_hist_2 <-  function(param_df,true_param,nsim){
  library(tidyverse)
  
   data_long <- param_df %>%
    pivot_longer(cols = c(rate, shape, beta_educ, beta_sex), 
                 names_to = "parameter", 
                 values_to = "value")
  
  true_param_long <- true_param %>%
    pivot_longer(cols = c(rate, shape, beta_educ, beta_sex), 
                 names_to = "parameter", 
                 values_to = "true_value")
  
  # Merge true parameter values with simulation data
  data_long <- data_long %>%
    left_join(true_param_long, by = c("trans", "parameter"))
  
  parameters <- unique(data_long$parameter)
  plot_list <- list()
  
  for (param in parameters) {
    plot <- ggplot(filter(data_long, parameter == param), aes(x = value, fill = model)) +
      geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
      geom_vline(aes(xintercept = true_value), color = "red", linetype = "dashed", linewidth = 1) +
      facet_wrap(~ trans, scales = "free") +  # Allow free scales for each transition
      theme_bw() +
      labs(
        title = paste("Histogram of", param, "Grouped by Transition and Model"),
        x = paste(param, "Value"),
        y = "Count"
      ) +
      theme(
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
      ) +
      scale_fill_brewer(palette = "Set2")  # color palette
    
    plot_list[[param]] <- plot
  }
  return(plot_list)
}

dotplot_rel_bias <- function(df, true_param,nsim){
  data_long <- df %>%
    pivot_longer(cols = c(rate, shape, beta_educ,beta_sex), 
                 names_to = "parameter", 
                 values_to = "value")
  
  true_param_long <- true_param %>%
    pivot_longer(cols = c(rate, shape, beta_educ, beta_sex), 
                 names_to = "parameter", 
                 values_to = "true_value")
  
  data_long <- data_long %>%
    left_join(true_param_long, by = c("trans", "parameter"))
  
  data_summary <- data_long %>%
    group_by(model, trans, parameter) %>%
    summarise(
      n = n(), 
      rel_bias = mean((value - true_value)/abs(true_value), na.rm = TRUE),
      se_rel_bias = sd((value - true_value)/abs(true_value), na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  # dot plot of mean values with error bars
  ggplot(data_summary, aes(x = model, y = rel_bias, color = trans)) +
    geom_point(size = 3) +  
    geom_errorbar(aes(ymin = rel_bias - se_rel_bias, ymax =  rel_bias + se_rel_bias), width = 0.2) +  
    geom_hline(aes(yintercept = 0), color = "gray", linetype = "dashed", linewidth = 1) +  
    facet_wrap(~ parameter, scales = "free_y", nrow = 2) +
    theme_bw() +
    labs(
      title = paste("Relative Bias by models and transitions. N=", nsim),
      x = "Model Name",
      y = "Relative Bias"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_manual(values = c( "blue", "purple", "pink"))
  
}

dotplot_rel_bias_misc <- function(df, true_param,nsim, ST){
  study_type <- ifelse(ST=="PS","Population based study", "Irregular visits")
  data_long <- df %>%
    pivot_longer(cols = c(rate, shape, beta_educ,beta_sex), 
                 names_to = "parameter", 
                 values_to = "value")
  
  true_param_long <- true_param %>%
    pivot_longer(cols = c(rate, shape, beta_educ, beta_sex), 
                 names_to = "parameter", 
                 values_to = "true_value")
  
  data_long <- data_long %>%
    left_join(true_param_long, by = c("trans", "parameter"))
  
  data_summary <- data_long %>%
    group_by(model, trans, parameter) %>%
    summarise(
      n = n(), 
      rel_bias = mean((value - true_value)/abs(true_value), na.rm = TRUE),
      se_rel_bias = sd((value - true_value)/abs(true_value), na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(misc = ifelse(model %in% c("IV_ApproxTIHMM","PS_ApproxTIHMM", "IV_TIHMM", "PS_TIHMM"), "Misclassification Included", "Misclassification Excluded"))
  
  etes_means <- data_summary %>%
    filter(model == "ETES") %>%
    dplyr::select(trans, parameter, ETES_mean = rel_bias) 
  
  data_summary <- data_summary %>%
    left_join(etes_means, by = c("trans", "parameter"))
  
  data_summary <- data_summary %>%
    filter(model != "ETES")
  
  data_summary <- data_summary %>%
    filter(
      (ST== "PS" & grepl("^PS_", model)) |
        (ST == "IV" & grepl("^IV_", model)) |
        !(ST %in% c("PS", "IV")) # Keep all models for other study types
    )
  # dot plot of mean values with error bars
  ggplot(data_summary, aes(x = model, y = rel_bias, color = misc)) +
    geom_point(size = 3) +  
    geom_errorbar(aes(ymin = rel_bias - se_rel_bias, ymax =  rel_bias + se_rel_bias), width = 0.2) +  
    #geom_hline(aes(yintercept = 0), color = "gray", linetype = "dashed", linewidth = 1) +  
    geom_hline(aes(yintercept = ETES_mean, linetype = "Benchmark Avg Rel Bias"), color = "#6699CC", linewidth = 1) +  
    facet_wrap(~parameter + trans, scales = "free_y", nrow = 4) +
    #facet_grid(parameter ~ trans, scales = "free_y") +
    #facet_wrap(~ parameter, scales = "free_y", nrow = 2) +
    theme_bw() +
    labs(
      title = paste("Relative Bias by models and transitions. ",study_type," - Nsim = ", nsim),
      x = "Model Name",
      y = "Relative Bias",
      linetype = "Reference Line" 
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_manual(values = c(  "lightblue", "darkblue"))+
    scale_linetype_manual(values = c("Benchmark Avg Rel Bias" = "dashed"))
  
}


dotplot_rel_bias_1_2<- function(df, true_param,nsim,ST=NULL){
  df <- df %>% filter(trans=="1->2")
  data_long <- df %>%
    pivot_longer(cols = c(beta_educ,beta_sex), 
                 names_to = "parameter", 
                 values_to = "value")
  
  true_param_long <- true_param %>%
    pivot_longer(cols = c(beta_educ, beta_sex), 
                 names_to = "parameter", 
                 values_to = "true_value")
  
  data_long <- data_long %>%
    left_join(true_param_long, by = c("trans", "parameter"))
  

  data_summary <- data_long %>%
    group_by(model, trans, parameter) %>%
    summarise(
      n = n(), 
      rel_bias = mean((value - true_value)/abs(true_value), na.rm = TRUE),
      se_rel_bias = sd((value - true_value)/abs(true_value), na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(misc = ifelse(model %in% c("IV_ApproxTIHMM","PS_ApproxTIHMM", "IV_TIHMM", "PS_TIHMM"), "Misclassification Included", "Misclassification Excluded"))
  
  etes_means <- data_summary %>%
    filter(model == "ETES") %>%
    dplyr::select(trans, parameter, ETES_mean = rel_bias) 
  
  data_summary <- data_summary %>%
    left_join(etes_means, by = c("trans", "parameter"))
  
  data_summary <- data_summary %>%
    filter(model != "ETES")
  if (!is.null(ST)){
    data_summary <- data_summary %>%
      filter(model=="ETES" |
               (ST == "PS" & grepl("^PS_", model)) |
               (ST == "IV" & grepl("^IV_", model)) |
               !(ST %in% c("PS", "IV")) # Keep all models for other study types
      )
  }
  
  data_summary <- data_summary %>%
    mutate(model = sub("^[^_]+_", "", model))
  
   ggplot(data_summary, aes(x = model, y = rel_bias, color = misc)) +
    geom_point(size = 3) +  
    geom_errorbar(aes(ymin = rel_bias - se_rel_bias, ymax =  rel_bias + se_rel_bias), width = 0.2) +  
    #geom_hline(aes(yintercept = 0), color = "gray", linetype = "dashed", linewidth = 1) +  
    geom_hline(aes(yintercept = ETES_mean, linetype = "Benchmark Avg Rel Bias"), color = "#6699CC", linewidth = 1) +  
    facet_wrap(~parameter, scales = "free_y") +
    #facet_grid(parameter ~ trans, scales = "free_y") +
    #facet_wrap(~ parameter, scales = "free_y", nrow = 2) +
    theme_bw() +
    labs(
      title = paste("Relative Bias in transitions from Mild to Complex (1->2). N=", nsim),
      x = "Model Name",
      y = "Relative Bias",
      linetype = "Reference Line" 
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    scale_color_manual(values = c(  "lightblue", "darkblue"))+
    scale_linetype_manual(values = c("Benchmark Avg Rel Bias" = "dashed"))
  
}


plot_est_stability <- function(param_df, true_param ,model_name, nsim){
  
  data_long <- param_df %>%
    pivot_longer(cols = c(rate, shape, beta_educ,beta_sex), 
                 names_to = "parameter", 
                 values_to = "value")
  
  true_param_long <- true_param %>%
    pivot_longer(cols = c(rate, shape, beta_educ, beta_sex), 
                 names_to = "parameter", 
                 values_to = "true_value")
  
  data_long <- data_long %>%
    left_join(true_param_long, by = c("trans", "parameter"))
  
  ###
  stability_df <- data_long %>%
    filter(model == model_name) %>%
    group_by(parameter, trans) %>%
    arrange(dataset_id) %>%
    mutate(avg_estimate = cummean(value),
           n_datasets = row_number())
  
  ggplot(stability_df, aes(x = n_datasets, y = avg_estimate, color = trans)) +
    geom_line() +
    geom_hline(aes(yintercept = true_value, color = trans), linetype = "dashed") +
    facet_wrap(~ parameter, scales = "free_y") +
    labs(title = paste("Stability of Estimates for", model_name, ". Nsim=",nsim),
         x = "Number of Datasets Included",
         y = "Average Estimated Value") +
    theme_minimal()
  
}
plot_est_stability2 <- function(param_df, true_param ,model_name, nsim){
  
  data_long <- param_df %>%
    pivot_longer(cols = c(rate, shape, beta_educ,beta_sex), 
                 names_to = "parameter", 
                 values_to = "value")
  
  true_param_long <- true_param %>%
    pivot_longer(cols = c(rate, shape, beta_educ, beta_sex), 
                 names_to = "parameter", 
                 values_to = "true_value")
  
  data_long <- data_long %>%
    left_join(true_param_long, by = c("trans", "parameter"))
  
  ###
  stability_df <- data_long %>%
    filter(model == model_name) %>%
    group_by(parameter, trans) %>%
    arrange(dataset_id) %>%
    mutate(avg_estimate = cummean(value),
           n_datasets = row_number())%>%
    filter(n_datasets >=25)
  
  ggplot(stability_df, aes(x = n_datasets, y = avg_estimate, color = trans)) +
    geom_line() +
    geom_hline(aes(yintercept = true_value, color = trans), linetype = "dashed") +
    facet_wrap(trans~ parameter, scales = "free_y") +
    labs(title = paste("Stability of Estimates for", model_name, ". Nsim=",nsim),
         x = "Number of Datasets Included",
         y = "Average Estimated Value") +
    scale_y_continuous(labels = function(x) signif(x, 3), 
                       expand = expansion(mult = 0.5)) +
  
    theme_minimal()
  
}

plot_cov_stability <- function(param_df, true_param ,model_name, nsim){
  
  data_long <- param_df %>%
    pivot_longer(cols = c(rate, shape, beta_educ,beta_sex), 
                 names_to = "parameter", 
                 values_to = "value")%>%
    pivot_longer(cols = c(rate_se, shape_se, beta_educ_se, beta_sex_se),
               names_to = "se_param", values_to = "se_value") %>%
    mutate(se_param = sub("_se$", "", se_param)) %>%
    filter(parameter == se_param) %>%
    dplyr::select(-se_param)
    
  true_param_long <- true_param %>%
    pivot_longer(cols = c(rate, shape, beta_educ, beta_sex), 
                 names_to = "parameter", 
                 values_to = "true_value")
  
  data_long <- data_long %>%
    left_join(true_param_long, by = c("trans", "parameter"))
  
  ###
  stability_df <- data_long %>%
    filter(model == model_name) %>%
    group_by(parameter, trans) %>%
    arrange(dataset_id) %>%
    mutate(coverage = sapply(1:n(), function(i) {
      calculate_coverage(value[1:i], se_value[1:i], true_value[1:i]) }),
           n_datasets = row_number()) %>%
    filter(n_datasets >=25)
  print(stability_df)
  ggplot(stability_df, aes(x = n_datasets, y = coverage, color = trans)) +
    geom_line() +
    #geom_hline(aes(yintercept = true_value, color = trans), linetype = "dashed") +
    facet_wrap(~ parameter, scales = "free_y") +
    labs(title = paste("Stability of Coverage for", model_name, ". Nsim=",nsim),
         x = "Number of Datasets Included",
         y = "Coverage") +
    theme_minimal()
  
}

calculate_coverage <- function(estimates, se, true_param, confidence_level = 0.95) {
  alpha <- 1 - confidence_level
  z <- qnorm(1 - alpha / 2)  
  
  lower_ci <- estimates - z * se
  upper_ci <- estimates + z * se
  within_ci <- (true_param >= lower_ci) & (true_param <= upper_ci)
  
  coverage <- mean(within_ci, na.rm = TRUE)
  
  
  return(coverage)
}
