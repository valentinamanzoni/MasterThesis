library(cowplot)  # For arranging plots
# obtain approximate gompertz parrameters from models fitted with msm
baseline_param <- function(model, model_name){
  rate<- numeric()
  rate_se <-numeric()
  shape <- numeric()
  shape_se <- numeric()
  trans <-numeric()
  haz <- hazards_snack(model,SE=TRUE, b.covariates = list(Age = 60, educ_el = 0, dm_sex = 0,no_pa=0,life_alone=0, heavy_alcool=0, if_ever_smoke=0,fin_strain_early=1,finstrain_dummy=0,sei_long_cat_dummy=0), no.years = 40)
  ntrans <-  sum(model$qmodel$imatrix)
  trans <- rename_trans(names(haz))
  min_age <- min(model$data[[1]]$Age)
  max_age <- max(model$data[[1]]$Age)
  
  # SE for rate and shape are computed with bootstrapping
  n_boot <- 1000
  for (i in 1:ntrans){
    age_grid <- seq(min_age, max_age , length.out = length(unlist(haz[[1]]$hazard)))
    y <- log(as.numeric(unlist(haz[[i]]$hazard)))
    reg_model <- lm(y ~ age_grid)
    rate <- c(rate,reg_model$coefficients[1])
    shape<- c(shape, reg_model$coefficients[2])
    haz_se <- unlist(haz[[i]]$hazardSE)
    
    bootstrap_results <- replicate(n_boot, {
      y_noisy <- log(as.numeric(unlist(haz[[i]]$hazard))) + rnorm(length(y), 0, haz_se)
      indices <- sample(1:length(y), replace = TRUE)
      y_boot <- y_noisy[indices]
      age_boot <- age_grid[indices]
      
      boot_model <- lm(y_boot ~ age_boot)
      boot_model$coefficients
    })
    
    rate_se <- c(rate_se, sd(bootstrap_results[1, ]))
    shape_se <- c(shape_se, sd(bootstrap_results[2, ]))
  }
  base_est <- data.frame(
    model = rep(model_name, length(trans)),
    trans = trans,
    rate = rate,
    rate_se = rate_se,
    shape = shape,
    shape_se = shape_se
  )
  return(base_est)
}

# plots: 

plot_cum_hazard_snack <- function(pm_df, t_vals, nsim=NULL) {
  
  for (i in unique(pm_df$trans)) {
    mp <- pm_df %>% filter(trans == i)
    
    model_ids <- unique(mp$model)
    model_hazards <- list()
    
    for (j in seq_along(model_ids)) {
      model_data <- mp %>% filter(model == model_ids[j])
      shape_param <- model_data$shape
      rate_param <- model_data$rate
      
      if (length(shape_param) == 1 && length(rate_param) == 1) {
        model_hazards[[j]] <- tibble(
          t = t_vals,
          hazard = Hgompertz(t_vals, shape_param, exp(rate_param), log = FALSE),
          model = paste0("Model ", model_ids[j])
        )
      }
    }
    
    
    hazard_data <-bind_rows(model_hazards)
    model_colors <- c(setNames(rainbow(length(model_ids)), paste0("Model ", model_ids)))
    model_linetypes <- c( setNames(rep("dashed", length(model_ids)), paste0("Model ", model_ids)))
    
    p <- ggplot(hazard_data, aes(x = t, y = hazard, color = model, linetype = model)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = model_colors) +
      scale_linetype_manual(values = model_linetypes) +
      labs(
        title = paste0("Cumulative Gompertz Hazard Function - Transition ", i),
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
plot_hazard_snack<-function(pm_df,  t_vals, nsim=NULL, ST=NULL) {
  
  plot_list <- list()  # Store individual plots
  
  for (i in unique(pm_df$trans)) {
    mp <- pm_df %>% filter(trans == i)
    
    model_ids <- unique(mp$model)
    
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
      shape_param <- model_data$shape
      rate_param <- model_data$rate
      
      if (length(shape_param) == 1 && length(rate_param) == 1) {
        model_hazards[[j]] <- tibble(
          t = t_vals,
          hazard = hgompertz(t_vals, shape_param, exp(rate_param), log = FALSE),
          model = paste0("Model ", model_ids[j])
        )
      }
    }
    
    hazard_data <- bind_rows(bind_rows(model_hazards))
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
    draw_label(paste("Gompertz Hazard Function Comparison"), fontface = "bold", size = 16, x = 0.5, y = 0.95) + 
    draw_plot(combined_plot, x = 0, y = 0, width = 1, height = 0.9)
  
  print(final_plot)
  
}
plot_hazard_snack_ci <- function(pm_df, t_vals, nsim = NULL, ST = NULL) {
  plot_list <- list()
  
  for (i in unique(pm_df$trans)) {
    mp <- pm_df %>% filter(trans == i)
    model_ids <- unique(mp$model)
    
    if (!is.null(ST)) {
      mp <- mp %>%
        filter(
          model == "ETES" |
            (ST == "PS" & grepl("^PS_", model)) |
            (ST == "IV" & grepl("^IV_", model)) |
            !(ST %in% c("PS", "IV"))
        )
    }
    
    model_ids <- unique(mp$model)
    model_hazards <- list()
    
    for (j in seq_along(model_ids)) {
      model_data <- mp %>% filter(model == model_ids[j])
      shape_param <- model_data$shape
      rate_param <- model_data$rate
      shape_se <- model_data$shape_se
      rate_se <- model_data$rate_se
      
      if (length(shape_param) == 1 && length(rate_param) == 1) {
        # Hazard and Delta Method for CIs
        hazard_vals <- hgompertz(t_vals, shape_param, exp(rate_param), log = FALSE)
        hazard_se <- hazard_vals * sqrt(rate_se^2 + (t_vals^2 * shape_se^2))
        
        model_hazards[[j]] <- tibble(
          t = t_vals,
          hazard = hazard_vals,
          lower = hazard_vals - 1.96 * hazard_se,
          upper = hazard_vals + 1.96 * hazard_se,
          model = paste0("Model ", model_ids[j])
        )
      }
    }
    
    hazard_data <- bind_rows(model_hazards)
    model_colors <- c("True" = "black", setNames(rainbow(length(model_ids)), paste0("Model ", model_ids)))
    model_linetypes <- c("True" = "solid", setNames(rep("dashed", length(model_ids)), paste0("Model ", model_ids)))
    
    p <- ggplot(hazard_data, aes(x = t, y = hazard, color = model, linetype = model)) +
      # Plot CI as ribbon
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = 0.2, color = NA) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = model_colors) +
      scale_fill_manual(values = model_colors) +
      scale_linetype_manual(values = model_linetypes) +
      labs(
        title = paste0("Gompertz Hazard Function - Transition ", i),
        x = "Age",
        y = "Hazard Function",
        color = "Model",
        linetype = "Model",
        fill = "Model"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(size = 12, hjust = 0.5)
      )
    
    plot_list[[as.character(i)]] <- p
  }
  
  legend <- get_legend(plot_list[[1]] + theme(legend.position = "right"))
  plots_no_legend <- lapply(plot_list, function(p) p + theme(legend.position = "none"))
  
  combined_plot <- plot_grid(
    plots_no_legend[[1]], plots_no_legend[[2]],
    plots_no_legend[[3]], legend,
    nrow = 2, ncol = 2,
    rel_widths = c(1, 1), rel_heights = c(1, 1)
  )
  
  final_plot <- ggdraw() +
    draw_label("Gompertz Hazard Function Comparison", fontface = "bold", size = 16, x = 0.5, y = 0.95) +
    draw_plot(combined_plot, x = 0, y = 0, width = 1, height = 0.9)
  
  print(final_plot)
}

plot.nhm.mine <- function (x, what = "probabilities", time0 = 0, state0 = 1, times = NULL, 
                           covvalue = NULL, ci = TRUE, sim = FALSE, coverage = 0.95, 
                           B = 1000, rtol = 1e-06, atol = 1e-06, main_arg = NULL, xlab = "Time") {
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  model <- x
  qmatrix.nhm <- function (object, time0 = 0, times = NULL, covvalue = NULL, ci = TRUE, 
                           sim = FALSE, coverage = 0.95, B = 1000) {
    model <- object
    if (is.null(times)) {
      times <- seq(time0, model$maxtime, length.out = 100)
    }
    npar <- model$model_object$nparQ
    par <- model$par[1:npar]
    ncov <- model$ncov
    fisher <- model$hess
    intens <- model$intens
    nstate <- model$nstate
    trans <- model$model_object$trans
    R <- nstate
    n1 <- (array(paste(rep(1:R, R), rep(1:R, each = R), sep = "->"), 
                 c(R, R)))
    n1 <- n1[c(trans) != 0]
    if (is.null(fisher)) 
      fisher <- model$fisher
    if (is.null(fisher) & ci) 
      stop("Fisher information required for confidence intervals")
    if (ci & model$singular) {
      warning("Cannot compute confidence intervals with a singular Hessian.")
      ci <- FALSE
    }
    if (is.null(covvalue) & ncov > 0) {
      covvalue <- model$covmeans
    }
    if (length(covvalue) != model$ncov) 
      stop("covvalue must be a vector of length equal to the number of covariates")
    nstate < model$nstate
    Q <- sapply(times, function(t) intens(t, covvalue, par))
    q1 <- t(sapply(Q[1, ], function(x) x))
    intensities <- data.frame(q1[, c(trans) != 0])
    names(intensities) <- n1
    row.names(intensities) <- round(times, 3)
    if (ci) {
      if (sim) {
        sig <- solve(fisher)[1:npar, 1:npar]
        sig <- 0.5 * (sig + t(sig))
        newpar <- rmvnorm(n = B, mean = par, sigma = sig)
        intB <- array(0, c(length(times), length(n1), B))
        for (i in 1:B) {
          Q <- sapply(times, function(t) intens(t, covvalue, 
                                                newpar[i, ]))
          q1 <- t(sapply(Q[1, ], function(x) x))
          intB[, , i] <- q1[, c(trans) != 0]
        }
        const <- 0.5 * (1 - coverage)
        intensL <- data.frame(apply(intB, c(1, 2), function(x) sort(x)[B * 
                                                                         const]))
        intensU <- data.frame(apply(intB, c(1, 2), function(x) sort(x)[B * 
                                                                         (1 - const)]))
      }
      else {
        covmat <- solve(fisher)[1:npar, 1:npar]
        dq1 <- array(sapply(Q[2, ], function(x) x), c(nstate^2, 
                                                      npar, length(times)))
        dq1 <- dq1[c(trans) != 0, , ]
        vars <- array(0, c(length(times), length(n1)))
        for (k in 1:length(times)) {
          for (l in 1:length(n1)) {
            vars[k, l] <- t(dq1[l, , k]) %*% covmat %*% 
              (dq1[l, , k])
          }
        }
        vars <- vars/(intensities^2 + 1 * (intensities == 
                                             0))
        const <- qnorm(1 - 0.5 * (1 - coverage))
        intensL <- intensities * exp(-const * vars^0.5)
        intensU <- intensities * exp(const * vars^0.5)
      }
      names(intensL) <- names(intensU) <- n1
      row.names(intensL) <- row.names(intensL) <- round(times, 
                                                        3)
    }
    else {
      intensL <- intensU <- NULL
    }
    list(times = times, intensities = intensities, lower = intensL, 
         upper = intensU)
  }
  if (ci & model$singular) {
    warning("Confidence intervals cannot be calculated due to a singular Hessian")
    ci <- FALSE
  }
  if (what == "probabilities") {
    if (is.null(main_arg)) 
      main_arg <- "Probability in state "
    preds <- predict.nhm.mine(model, time0, state0, times, covvalue, 
                              ci, sim, coverage, B, rtol, atol)
    nstate <- model$nstate
    prevalence <- preds$probabilities
    times <- preds$times
    prevL <- preds$lower
    prevU <- preds$upper
    par(mfrow = c(2, ceiling(nstate/2)))
    for (i in 1:nstate) {
      plot(times, prevalence[, i], type = "l", xlab = xlab, 
           ylab = "Probability", main = paste(main_arg, 
                                              i, sep = ""), ylim = c(0, 1))
      if (ci) {
        lines(times, prevL[, i], lty = 2)
        lines(times, prevU[, i], lty = 2)
      }
      else {
        prevL <- prevU <- NULL
      }
    }
  }
  else {
    if (what == "intensities") {
      if (is.null(main_arg)) 
        main_arg <- "Intensity "
      preds <- qmatrix.nhm(model, time0, times, covvalue, 
                           ci, sim, coverage, B)
      ntrans <- dim(preds$intensities)[2]
      intensnames <- names(preds$intensities)
      intensities <- preds$intensities
      times <- preds$times
      prevL <- preds$lower
      prevU <- preds$upper
      if (ci) {
        maxintens <- apply(prevU, 2, stats::quantile, 
                           0.98)
      }
      else {
        maxintens <- apply(intensities, 2, stats::quantile, 
                           0.98)
      }
      #par(mfrow = c(2, ceiling(ntrans/2)))
      plots <- list()
      for (i in 1:ntrans) {
        plot_df <- data.frame(
          time = times,
          intensity = intensities[, i],
          lower = if (ci) prevL[, i] else NA,
          upper = if (ci) prevU[, i] else NA,
          transition = intensnames[i]
        )
        plot_df$type <- ifelse(is.na(plot_df$lower), "Hazard Function", "95% CI")
        
        p <- ggplot(plot_df, aes(x = time, y = intensity)) +
          geom_line(aes(color = "Hazard Function"), linewidth = 0.5) +
          labs(x = xlab, y = "Hazard", title = paste(main_arg, intensnames[i])) +
          scale_color_manual(name = "Legend", values = c("Hazard Function" = "blue", "95% CI" = "lightblue")) +
          scale_fill_manual(name ="Region", values = c("Hazard Function" = "blue", "95% CI" = "lightblue")) +
          theme_minimal() +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "lightgray", fill = NA, linewidth = 0.5),
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
            legend.position = "bottom" 
          )
        
        
        
        # Add ribbon only if CI is available
        if (ci) {
          p <- p + 
            geom_ribbon(aes(ymin = lower, ymax = upper, fill = type), alpha = 0.4) 
        }
        
        plots[[i]] <- p 
      }
      legend <- get_legend(plots[[1]] + theme(legend.position = "right"))
      plots_no_legend <- lapply(plots, function(p) p + theme(legend.position = "none"))
      
      combined_plot <- plot_grid(
        plots_no_legend[[1]], plots_no_legend[[2]], 
        plots_no_legend[[3]], legend, 
        nrow = 2, ncol = 2, 
        rel_widths = c(1, 1), rel_heights = c(1, 1)
      )
      final_plot <- ggdraw() +
        draw_plot(combined_plot, x = 0, y = 0, width = 1, height = 0.9)
      
      print(final_plot)
    }
    if (what != "intensities") 
      stop("Only probabilities or intensities can be plotted.")
  }
}

predict.nhm.mine <- function (object, time0 = 0, state0 = 1, times = NULL, covvalue = NULL, 
                              ci = TRUE, sim = FALSE, coverage = 0.95, B = 1000, rtol = 1e-06, 
                              atol = 1e-06) {
  model <- object
  if (time0 == 0 & model$model_object$type == "weibull") 
    time0 <- 1e-05
  if (is.null(times)) {
    times <- seq(time0, model$maxtime, length.out = 100)
  }
  genintens.nhm <- function (t, y, parms, intens) {
    npar <- parms[1]
    ncov <- parms[2]
    nstate <- parms[3]
    p <- array(y, dim = c(nstate, nstate))
    if (ncov > 0) {
      z <- parms[4:(3 + ncov)]
    }
    else {
      z <- NULL
    }
    x <- parms[(4 + ncov):(3 + ncov + npar)]
    Q <- intens(t, z = z, x = x)
    q <- Q$q
    vec <- c((p %*% q))
    list(vec, NULL)
  }
  genintens_deriv.nhm <- function (t, y, parms, intens) {
    npar <- parms[1]
    ncov <- parms[2]
    nstate <- parms[3]
    p <- array(y[1:(nstate^2)], dim = c(nstate, nstate))
    pp <- array(0, c(nstate, nstate, npar))
    for (j in 1:npar) {
      pp[, , j] <- array(y[(1 + j * nstate^2):((j + 1) * nstate^2)], 
                         dim = c(nstate, nstate))
    }
    if (ncov > 0) {
      z <- parms[4:(3 + ncov)]
    }
    else {
      z <- NULL
    }
    x <- parms[(4 + ncov):(3 + ncov + npar)]
    Q <- intens(t, z = z, x = x)
    q <- Q$q
    qp <- Q$qp
    vec <- c((p %*% q))
    for (j in 1:npar) {
      vec <- c(vec, c(p %*% qp[, , j] + pp[, , j] %*% q))
    }
    list(vec, NULL)
  }
  npar <- model$model_object$nparQ
  par <- model$par[1:npar]
  ncov <- model$ncov
  fisher <- model$hess
  intens <- model$intens
  nstate <- model$nstate
  if (is.null(fisher)) 
    fisher <- model$fisher
  if (is.null(fisher) & ci) 
    stop("Fisher information required for confidence intervals")
  if (model$singular & ci) {
    warning("Cannot calculate confidence intervals due to singular Hessian")
    ci <- FALSE
  }
  if (is.null(covvalue) & ncov > 0) {
    covvalue <- model$covmeans
  }
  if (length(covvalue) != model$ncov) 
    stop("covvalue must be a vector of length equal to the number of covariates")
  nstate < model$nstate
  tcrit <- model$tcrit
  if (is.null(tcrit)) 
    tcrit <- max(times) + 2
  if (ncov > 0) {
    parms <- c(npar, ncov, nstate, covvalue, par)
  }
  else {
    parms <- c(npar, ncov, nstate, par)
  }
  if (min(times) != time0) 
    times <- c(time0, times)
  res <- deSolve::lsoda(y = c(diag(nstate)), times = times, 
                        func = genintens.nhm, tcrit = tcrit, rtol = rtol, atol = atol, 
                        parms = parms, intens = intens)
  res2 <- array(res[, 2:(nstate^2 + 1)], c(length(times), nstate, 
                                           nstate))
  prevalence <- res2[, state0, ]
  if (ci) {
    if (sim) {
      sig <- solve(fisher)[1:npar, 1:npar]
      sig <- 0.5 * (sig + t(sig))
      newpar <- rmvnorm(n = B, mean = par, sigma = sig)
      prevB <- array(0, c(length(times), nstate, B))
      for (i in 1:B) {
        if (ncov > 0) {
          parms <- c(npar, ncov, nstate, covvalue, newpar[i, 
          ])
        }
        else {
          parms <- c(npar, ncov, nstate, newpar[i, ])
        }
        res <- deSolve::lsoda(y = c(diag(nstate)), times = times, 
                              func = genintens.nhm, tcrit = tcrit, rtol = rtol, 
                              atol = atol, parms = parms, intens = intens)
        res2 <- array(res[, 2:(nstate^2 + 1)], c(length(times), 
                                                 nstate, nstate))
        prevB[, , i] <- res2[, state0, ]
      }
      const <- 0.5 * (1 - coverage)
      prevL <- apply(prevB, c(1, 2), function(x) sort(x)[B * 
                                                           const])
      prevU <- apply(prevB, c(1, 2), function(x) sort(x)[B * 
                                                           (1 - const)])
    }
    else {
      covmat <- solve(fisher)[1:npar, 1:npar]
      res <- deSolve::lsoda(y = c(c(diag(nstate)), rep(0, 
                                                       nstate^2 * npar)), times = times, func = genintens_deriv.nhm, 
                            tcrit = tcrit, rtol = rtol, atol = atol, parms = parms, 
                            intens = intens)
      p0 <- array(res[, 2:(nstate^2 + 1)], c(length(times), 
                                             nstate, nstate))
      lik <- 0
      dp0 <- array(0, c(length(times), nstate, nstate, 
                        npar))
      for (l in 1:npar) {
        q0 <- array(res[, (2 + l * nstate^2):((l + 1) * 
                                                nstate^2 + 1)], c(length(times), nstate, nstate))
        dp0[, , , l] <- q0
      }
      dprev <- dp0[, state0, , ]
      vars <- array(0, c(length(times), nstate))
      for (k in 1:length(times)) {
        for (l in 1:nstate) {
          vars[k, l] <- t(dprev[k, l, ]) %*% covmat %*% 
            (dprev[k, l, ])
        }
      }
      vars <- vars/(prevalence * (1 - prevalence) + 1 * 
                      (prevalence == 0) + 1 * (prevalence == 1))^2
      lprevalence <- log(prevalence) - log(1 - prevalence)
      const <- qnorm(1 - 0.5 * (1 - coverage))
      prevL <- 1/(1 + exp(-(lprevalence - const * vars^0.5)))
      prevU <- 1/(1 + exp(-(lprevalence + const * vars^0.5)))
    }
  }
  else {
    prevL <- prevU <- NULL
  }
  list(times = times, probabilities = prevalence, lower = prevL, 
       upper = prevU)
}

