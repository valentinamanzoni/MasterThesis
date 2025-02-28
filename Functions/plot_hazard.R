
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)  # For arranging plots

plot_cum_hazard_gg <- function(pm_df, true_param, t_vals, nsim=NULL) {
  
  for (i in unique(true_param$trans)) {
    tp <- true_param %>% filter(trans == i)
    mp <- pm_df %>% filter(trans == i, parameter %in% c("rate", "shape"))
    
    true_hazard_vals <- tibble(
      t = t_vals,
      hazard = Hgompertz(t_vals, tp$shape, exp(tp$rate), log = FALSE),
      model = "True"
    )
    

    model_ids <- unique(mp$model)
    model_hazards <- list()
    
    for (j in seq_along(model_ids)) {
      model_data <- mp %>% filter(model == model_ids[j])
      shape_param <- model_data %>% filter(parameter == "shape") %>% pull(estimate)
      rate_param <- model_data %>% filter(parameter == "rate") %>% pull(estimate)
      
      if (length(shape_param) == 1 && length(rate_param) == 1) {
        model_hazards[[j]] <- tibble(
          t = t_vals,
          hazard = Hgompertz(t_vals, shape_param, exp(rate_param), log = FALSE),
          model = paste0("Model ", model_ids[j])
        )
      }
    }
    

    hazard_data <- bind_rows(true_hazard_vals, bind_rows(model_hazards))
    model_colors <- c("True" = "black", setNames(rainbow(length(model_ids)), paste0("Model ", model_ids)))
    model_linetypes <- c("True" = "solid", setNames(rep("dashed", length(model_ids)), paste0("Model ", model_ids)))
    
    p <- ggplot(hazard_data, aes(x = t, y = hazard, color = model, linetype = model)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = model_colors) +
      scale_linetype_manual(values = model_linetypes) +
      labs(
        title = paste0("Cumulative Gompertz Hazard Function - Transition ", i, " - Nsim:", nsim),
        x = "Time",
        y = "Cumulative Hazard Function",
        color = "Model",
        linetype = "Model"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top")
    
    print(p)  
  }
}

plot_hazard_gg <- function(pm_df, true_param, t_vals, nsim=NULL, ST=NULL) {
  
  for (i in unique(true_param$trans)) {
    tp <- true_param %>% filter(trans == i)
    mp <- pm_df %>% filter(trans == i, parameter %in% c("rate", "shape"))
    
    true_hazard_vals <- tibble(
      t = t_vals,
      hazard = hgompertz(t_vals, tp$shape, exp(tp$rate), log = FALSE),
      model = "True"
    )
    
    if (!is.null(ST)){
      mp <- mp %>%
        filter(
          (ST== "PS" & grepl("^PS_", model)) |
            (ST == "IV" & grepl("^IV_", model)) |
            grepl("ETES", model) | 
            !(ST %in% c("PS", "IV")) # Keep all models for other study types
        )
    }

    
    model_ids <- unique(mp$model)
    model_hazards <- list()
    
    
    for (j in seq_along(model_ids)) {
      model_data <- mp %>% filter(model == model_ids[j])
      shape_param <- model_data %>% filter(parameter == "shape") %>% pull(estimate)
      rate_param <- model_data %>% filter(parameter == "rate") %>% pull(estimate)
      
      if (length(shape_param) == 1 && length(rate_param) == 1) {
        model_hazards[[j]] <- tibble(
          t = t_vals,
          hazard = hgompertz(t_vals, shape_param, exp(rate_param), log = FALSE),
          model = paste0("Model ", model_ids[j])
        )
      }
    }
    
    
    hazard_data <- bind_rows(true_hazard_vals, bind_rows(model_hazards))
    model_colors <- c("True" = "black", setNames(rainbow(length(model_ids)), paste0("Model ", model_ids)))
    model_linetypes <- c("True" = "solid", setNames(rep("dashed", length(model_ids)), paste0("Model ", model_ids)))
    
    p <- ggplot(hazard_data, aes(x = t, y = hazard, color = model, linetype = model)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = model_colors) +
      scale_linetype_manual(values = model_linetypes) +
      labs(
        title = paste0("Gompertz Hazard Function - Transition ", i, " - Nsim:", nsim),
        x = "Age",
        y = "Hazard Function",
        color = "Model",
        linetype = "Model"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top")
    
    print(p)  
  }
}

plot_hazard_gg_2<-function(pm_df, true_param, t_vals, nsim=NULL, ST=NULL) {

  
  plot_list <- list()  # Store individual plots
  
  for (i in unique(true_param$trans)) {
    tp <- true_param %>% filter(trans == i)
    mp <- pm_df %>% filter(trans == i, parameter %in% c("rate", "shape"))
    
    true_hazard_vals <- tibble(
      t = t_vals,
      hazard = hgompertz(t_vals, tp$shape, exp(tp$rate), log = FALSE),
      model = "True"
    )
    
    if (!is.null(ST)){
      mp <- mp %>%
        filter(model=="ETES" |
          (ST == "PS" & grepl("^PS_", model)) |
            (ST == "IV" & grepl("^IV_", model)) |
            !(ST %in% c("PS", "IV")) # Keep all models for other study types
        )
    }
    
    model_ids <- unique(mp$model)
    model_hazards <- list()
    
    for (j in seq_along(model_ids)) {
      model_data <- mp %>% filter(model == model_ids[j])
      shape_param <- model_data %>% filter(parameter == "shape") %>% pull(estimate)
      rate_param <- model_data %>% filter(parameter == "rate") %>% pull(estimate)
      
      if (length(shape_param) == 1 && length(rate_param) == 1) {
        model_hazards[[j]] <- tibble(
          t = t_vals,
          hazard = hgompertz(t_vals, shape_param, exp(rate_param), log = FALSE),
          model = paste0("Model ", model_ids[j])
        )
      }
    }
    
    hazard_data <- bind_rows(true_hazard_vals, bind_rows(model_hazards))
    model_colors <- c("True" = "black", setNames(rainbow(length(model_ids)), paste0("Model ", model_ids)))
    model_linetypes <- c("True" = "solid", setNames(rep("dashed", length(model_ids)), paste0("Model ", model_ids)))
    
    p <- ggplot(hazard_data, aes(x = t, y = hazard, color = model, linetype = model)) +
      geom_line(data = hazard_data %>% filter(model == "True"), 
                aes(color = model, linetype = model), linewidth = 1.0) +
      # Plot other models on top
      geom_line(data = hazard_data %>% filter(model != "True"), 
                aes(color = model, linetype = model), linewidth = 1) +
      scale_color_manual(values = model_colors) +
      scale_linetype_manual(values = model_linetypes) +
      labs(
        title = paste0("Gompertz Hazard Function - Transition ", i),
        x = "Age",
        y = "Hazard Function",
        color = "Model",
        linetype = "Model"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom",
            plot.title = element_text(size = 12,  hjust = 0.5)) # Set legend to bottom for easy extraction
    
    plot_list[[as.character(i)]] <- p  # Store plot in the list
  }
  
  # Extract legend from one plot
  legend <- get_legend(plot_list[[1]] + theme(legend.position = "right"))
  
  # Remove legends from individual plots
  plots_no_legend <- lapply(plot_list, function(p) p + theme(legend.position = "none"))
  
  combined_plot <- plot_grid(
    plots_no_legend[[1]], plots_no_legend[[2]], 
    plots_no_legend[[3]], legend, 
    nrow = 2, ncol = 2, 
    rel_widths = c(1, 1), rel_heights = c(1, 1)
  )
  final_plot <- ggdraw() +
    draw_label(paste("Gompertz Hazard Function Comparison", " - Nsim:", nsim), fontface = "bold", size = 16, x = 0.5, y = 0.95) + 
    draw_plot(combined_plot, x = 0, y = 0, width = 1, height = 0.9)
  
  print(final_plot)
  
}
plot_cum_hazard<- function(pm_df, true_param, t_vals, nsim=NULL){
  
  for (i in true_param$trans){
    tp <- true_param%>% filter(trans==i)
    mp <- pm_df%>% filter(trans==i,parameter %in% c("rate","shape"))
    hazard_vals_t <- Hgompertz(t_vals, tp$shape, exp(tp$rate) , log = FALSE)
    plot(t_vals, hazard_vals_t, type = "l", col = "black", lwd = 2, 
         ylab = "Cumulative Hazard Function", xlab = "x", main = paste0("Cumulative Gompertz Hazard Function - Transition:",i, " - Nsim:", nsim))
    colors <- rainbow(nrow(mp) / 2)
    model_ids <- unique(mp$model)  
    for (j in seq_along(model_ids)) {
      model_data <- mp %>% filter(model == model_ids[j])
      shape_param <- model_data %>% filter(parameter == "shape") %>% pull(estimate)
      rate_param <- model_data %>% filter(parameter == "rate") %>% pull(estimate)
      
      if (length(shape_param) == 1 && length(rate_param) == 1) {
        hazard_vals_m <- Hgompertz(t_vals, shape_param, exp(rate_param), log = FALSE)
        lines(t_vals, hazard_vals_m, col = colors[j], lwd = 2, lty = 2)  
      }
    }
    
    legend("topleft", legend = c("True", paste0("Model ", model_ids)), 
           col = c("black", colors), lty = c(1, rep(2, length(model_ids))), lwd = 2)
  }
}

eval_cum_hazard <- function(pm_df, true_param, t_vals, nsim=NULL, ST=NULL){
  N <- length(t_vals)
  eval_cm <- data.frame()
  for (i in unique(true_param$trans)){
    tp <- true_param %>% filter(trans == i)
    mp <- pm_df %>% filter(trans == i, parameter %in% c("rate", "shape"))
    
    true_hazard_vals <- tibble(
      t = t_vals,
      hazard = hgompertz(t_vals, tp$shape, exp(tp$rate), log = FALSE),
      model = "True"
    )
    
    if (!is.null(ST)){
      mp <- mp %>%
        filter(
          model=="ETES" |
          (ST== "PS" & grepl("^PS_", model)) |
            (ST == "IV" & grepl("^IV_", model)) |
            !(ST %in% c("PS", "IV")) # Keep all models for other study types
        )
    }
    
    
    model_ids <- unique(mp$model)
    
    
    for (j in seq_along(model_ids)) {
      model_data <- mp %>% filter(model == model_ids[j])
      shape_param <- model_data %>% filter(parameter == "shape") %>% pull(estimate)
      rate_param <- model_data %>% filter(parameter == "rate") %>% pull(estimate)
      
      if (length(shape_param) == 1 && length(rate_param) == 1) {
        model_hazard_vals <- tibble(
          t = t_vals,
          hazard = hgompertz(t_vals, shape_param, exp(rate_param), log = FALSE),
          model = paste0("Model ", model_ids[j])
        )
        # mean squared error
        mse <- mean((true_hazard_vals$hazard - model_hazard_vals$hazard)^2)
        # integrated absolute error
        iae <- sum(abs(true_hazard_vals$hazard - model_hazard_vals$hazard)) * diff(t_vals[1:2])
        ext_KL <- sum(model_hazard_vals$hazard*(true_hazard_vals$hazard/model_hazard_vals$hazard - log(true_hazard_vals$hazard/model_hazard_vals$hazard)-1) * diff(t_vals[1:2]))
        
        eval_cm <-  bind_rows(eval_cm, tibble(
          transition = i,
          model = model_ids[j],
          MSE = mse,
          IAE = iae,
          EXT_KL= ext_KL
        ))
      }
    }
  }
    return(eval_cm)
}

  
