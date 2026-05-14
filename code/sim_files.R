library(survival)
library(ggplot2)
library(patchwork)
library(dplyr)

sim_coverage <- function( N = 1000, 
                          reps = 1000, 
                          HR = 2, 
                          t_grid=c(seq(0.00001, 0.1, by = 0.00001), seq(0.1, 1, by = 0.001)), 
                          newpatient = data.frame(X1 = 0, X2 = 0) ){
  # Parameters
  set.seed(123)
  beta <- log(HR)
  true_S_pt1 <- exp(-exp(beta * newpatient$X1 + beta * newpatient$X2) * t_grid)
  
  # Storage
  coverage_matrix <- matrix(NA, nrow = reps, ncol = length(t_grid))
  neff_matrix     <- matrix(NA, nrow = reps, ncol = length(t_grid))
  
  for (m in 1:reps) {
    # Data Gen
    X1 <- rbinom(N, 1, 0.5)
    X2 <- rnorm(N, 0, 1)
    hazards <- exp(beta * X1 + beta * X2)
    TT <- rexp(N, rate = hazards)
    df <- data.frame(TT = TT, status = 1, X1 = X1, X2 = X2)
    
    # Fit & Predict
    fit <- coxph(Surv(TT, status) ~ X1 + X2, data = df)
    s_obj <- summary(survfit(fit, newdata = newpatient, conf.type = "log-log"), times = t_grid)
    
    # Coverage
    coverage_matrix[m, ] <- (true_S_pt1 <= s_obj$upper) & (true_S_pt1 >= s_obj$lower)
    
    # n.eff Calculation + NaN fix
    n_eff <- (1 - s_obj$surv) / (s_obj$surv * (s_obj$std.err^2))
    
    # Handling NaNs/Infs (Standard error is 0 before first event)
    first_valid <- which(!is.na(n_eff) & !is.infinite(n_eff))[1]
    if (!is.na(first_valid) && first_valid > 1) {
      n_eff[1:(first_valid - 1)] <- n_eff[first_valid]
    }
    neff_matrix[m, ] <- n_eff
  }
  
  # Prepare Data for ggplot
  plot_df <- data.frame(
    time = t_grid,
    coverage = colMeans(coverage_matrix, na.rm = TRUE),
    n_eff_mean = colMeans(neff_matrix, na.rm = TRUE),
    true = true_S_pt1
  )
  
  # --- GGPLOT CONSTRUCTION ---
  
  # Top Plot: Coverage
  p1 <- ggplot(plot_df, aes(x = time, y = coverage)) +
    geom_line(color = "steelblue", size = 1) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    scale_y_continuous(limits = c(0.85, 1.0)) +
    labs(title = paste0("Simulation Results for patient (X1=",newpatient$X1, ", X2=", newpatient$X2, ")"), 
         x = "Time",
         y = "95% CI Coverage") +
    theme_minimal()
  
  # Bottom Plot: Effective Sample Size
  p2 <- ggplot(plot_df, aes(x = time, y = n_eff_mean)) +
    geom_line(color = "darkgreen", size = 1) +
    labs(x = "Time", y = "Mean Effective Sample Size (n.eff)") +
    geom_hline(yintercept = 1000, linetype = "dashed", color = "red") +
    scale_y_continuous(limits = c(0, 2500)) +
    theme_minimal()
  
  p3 <- ggplot(plot_df, aes(x=time, y=true)) +
    geom_line(color = "purple", size = 1) +
    labs(x = "Time", y = "True Survival Probability") +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal()
  
  # Concatenate top-bottom
  final_plot <- p1 / p2 / p3
  return(final_plot)
}

sim_se_comparison <- function(N = 1000, 
                              reps = 1000, 
                              HR = 2, 
                              t_grid = c(seq(0.00001, 0.1, by = 0.00001), seq(0.1, 1, by = 0.001)), 
                              newpatient = data.frame(X1 = 0, X2 = 0)) {
  set.seed(123)
  beta <- log(HR)
  
  # Storage for survival estimates and their model-reported SEs
  surv_estimates <- matrix(NA, nrow = reps, ncol = length(t_grid))
  model_se      <- matrix(NA, nrow = reps, ncol = length(t_grid))
  
  for (m in 1:reps) {
    # Data Generation (Exponential Baseline Hazard = 1)
    X1 <- rbinom(N, 1, 0.5)
    X2 <- rnorm(N, 0, 1)
    hazards <- exp(beta * X1 + beta * X2)
    TT <- rexp(N, rate = hazards)
    df <- data.frame(TT = TT, status = 1, X1 = X1, X2 = X2)
    
    # Fit Cox Model
    fit <- coxph(Surv(TT, status) ~ X1 + X2, data = df)
    
    # Extract survival and SE at t_grid
    # We use type="survival" to ensure we get the standard error of the estimate
    s_obj <- summary(survfit(fit, newdata = newpatient), times = t_grid)
    
    surv_estimates[m, ] <- s_obj$surv
    model_se[m, ]      <- s_obj$std.err
  }
  
  # Calculate Comparison Metrics
  # Empirical SE is the standard deviation of the point estimates across reps
  emp_se <- apply(surv_estimates, 2, sd, na.rm = TRUE)
  # Mean Model SE is the average of the SEs reported by the Cox model
  mean_model_se <- sqrt(colMeans(model_se^2, na.rm = TRUE))
  
  plot_df <- data.frame(
    time = t_grid,
    Emp_SE = emp_se,
    Model_SE = mean_model_se,
    Ratio = mean_model_se / emp_se
  )
  
  # Plot 1: Direct Comparison
  p1 <- ggplot(plot_df, aes(x = time)) +
    geom_line(aes(y = Emp_SE, color = "Empirical SE"), size = 1) +
    geom_line(aes(y = Model_SE, color = "Mean Model SE"), linetype = "dashed", size = 1) +
    scale_color_manual(values = c("Empirical SE" = "black", "Mean Model SE" = "red")) +
    labs(title = "Model-Based vs. Empirical Standard Errors",
         subtitle = paste0("Patient: X1=", newpatient$X1, ", X2=", newpatient$X2),
         y = "Standard Error", color = "Type") +
    theme_minimal()
  
  # Plot 2: Ratio (Should be close to 1.0)
  p2 <- ggplot(plot_df, aes(x = time, y = Ratio)) +
    geom_line() +
    geom_hline(yintercept = 1, linetype = "dotted", color = "blue") +
    scale_y_continuous(limits = c(0.8, 1.2)) +
    labs(y = "Ratio (Model SE / MC SE)", x = "Time") +
    theme_minimal()
  
  return(p1 / p2)
}

extract_sf <- function(sf_obj, tx, prevalence){
  df <- sf_to_df(sf_obj, time_grid)
  df$tx <- factor(tx)
  df$prevalence <- prevalence
  df
}

run_one <- function(rep_id, scale, shape, maxtime, hr, n, time_grid,
                    model_type = c("standard", "robust")){
  
  print(paste("Running replicate ", rep_id))
  
  model_type <- match.arg(model_type)
  
  out_list <- list()
  
  prev_list <- exp(-scale*10^shape)
  
  for(i in seq_along(scale)){
    
    prev_label <- prev_list[i]
    
    # simulate
    dat <- simsurv(
      dist = "weibull",
      gammas = shape,
      lambdas = scale[i],
      x = data.frame(tx = c(rep(0, n/2), rep(1, n/2))),
      betas = c("tx" = log(hr)),
      maxt = maxtime
    )
    
    dat$tx <- c(rep(0, n/2), rep(1, n/2))
    
    # choose model
    if(model_type == "standard"){
      fit <- coxph(Surv(eventtime, status) ~ tx, data = dat)
      model_label <- "Standard"
    } else {
      fit <- coxph(Surv(eventtime, status) ~ tx, data = dat, robust = TRUE)
      model_label <- "Robust"
    }
    
    # survfit_n
    sf_n0 <- survfit_n(survfit(fit, newdata = data.frame(tx = 0)), fit, coef=TRUE, chaz=TRUE)
    sf_n1 <- survfit_n(survfit(fit, newdata = data.frame(tx = 1)), fit, coef=TRUE, chaz=TRUE)
    
    df_i <- bind_rows(
      extract_sf(sf_n0, 0, prev_label),
      extract_sf(sf_n1, 1, prev_label)
    )
    
    df_i$model <- model_label
    
    out_list[[i]] <- df_i
  }
  
  df_rep <- bind_rows(out_list)
  df_rep$replicate <- rep_id
  
  return(df_rep)
}

lambda <- function( surv, shape ){
  -log(surv)/10^(shape)
}

run_one_generalized <- function(rep_id, scale, shape, maxtime, hr_vector, 
                                covariates, time_grid, 
                                model_type = c("standard", "robust")) {
  
  model_type <- match.arg(model_type)
  out_list <- list()
  betas <- log(hr_vector)
  prev_list <- exp(-scale * 10^shape)
  
  # Dynamic formula: Surv ~ var1 + var2 + ...
  cov_names <- names(covariates)
  formula_str <- as.formula(paste("Surv(eventtime, status) ~", paste(cov_names, collapse = " + ")))
  
  for(i in seq_along(scale)) {
    prev_label <- prev_list[i]
    
    # Simulate
    dat <- simsurv(
      dist = "weibull", gammas = shape, lambdas = scale[i],
      x = covariates, betas = betas, maxt = maxtime
    )
    dat <- cbind(dat, covariates)
    
    # Fit
    fit <- coxph(formula_str, data = dat, robust = (model_type == "robust"))
    
    # Generate predictions for every unique covariate profile
    unique_scenarios <- distinct(covariates)
    sf_list <- list()
    
    for(j in 1:nrow(unique_scenarios)) {
      row_data <- unique_scenarios[j, , drop = FALSE]
      
      # Get survfit and convert to df using your sf_to_df
      sf_raw <- survfit(fit, newdata = row_data)
      sf_n   <- survfit_n(sf_raw, fit, coef = TRUE, chaz = TRUE)
      
      # Use sf_to_df directly to avoid the old tx/prevalence hardcoding
      df_temp <- sf_to_df(sf_n, time_grid)
      
      # Attach prevalence and the specific covariate values for this profile
      df_temp$prevalence <- prev_label
      df_temp <- cbind(df_temp, row_data) 
      
      sf_list[[j]] <- df_temp
    }
    
    df_i <- bind_rows(sf_list)
    df_i$model <- ifelse(model_type == "robust", "Robust", "Standard")
    out_list[[i]] <- df_i
  }
  
  df_rep <- bind_rows(out_list)
  df_rep$replicate <- rep_id
  return(df_rep)
}

prepare_simulation_summary <- function(df_raw) {
  # 1. Identify all non-output columns as covariates
  # Based on your error, columns like 'lower', 'upper', etc., are present
  meta_cols <- c("time", "surv", "n.eff", "n.risk", "prevalence", "replicate", 
                 "model", "lower", "upper", "n.uncensor", "n.lower", "n.upper", 
                 "std.err", "std.chaz", "std.coef")
  
  cov_cols <- setdiff(names(df_raw), meta_cols)
  
  # 2. Average across replicates
  df_summary <- df_raw %>%
    group_by(across(all_of(c(cov_cols, "prevalence", "time")))) %>%
    summarise(
      n.eff = mean(n.eff, na.rm = TRUE),
      surv = mean(surv, na.rm = TRUE),
      n.risk = mean(n.risk, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 3. Calculate Profile ID and Events per group
  df_summary <- df_summary %>%
    group_by(across(all_of(c(cov_cols, "prevalence")))) %>%
    mutate(
      profile_id = cur_group_id(),
      # Number of events = initial risk - current risk
      events = max(n.risk, na.rm = TRUE) - n.risk,
      # Find the max events for this group to avoid division by zero
      max_ev = max(events, na.rm = TRUE)
    ) %>%
    mutate(
      # Use a basic ifelse to avoid the strict size-matching error of if_else
      events_rate = ifelse(max_ev == 0, 0, events / max_ev)
    ) %>%
    select(-max_ev) %>% # Clean up helper column
    ungroup()
  
  # 4. Extract Max Points for the frontier
  max_points <- df_summary %>%
    group_by(profile_id) %>%
    filter(time == max(time)) %>%
    ungroup()
  
  return(list(summary = df_summary, max_points = max_points, cov_cols = cov_cols))
}

plot_survival_simplified <- function(summary_list, y_limit = 1600, y_line = 1000, time_scale = "linear") {
  library(ggplot2)
  library(patchwork)
  
  df_summary <- summary_list$summary
  max_points <- summary_list$max_points
  
  # Define the simple gradient: Light to Dark
  # High survival (1.0) = Light Blue, Low survival (0.0) = Dark Blue
  simple_gradient <- scale_color_gradient(low = "#132B43", high = "#56B1F7", name = "Survival")
  
  # Base theme for consistency
  base_theme <- theme_minimal() +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  # --- Panel 1: Time ---
  p1 <- ggplot(df_summary, aes(x = time, y = n.eff, group = profile_id, color = surv)) +
    geom_line(linewidth = 0.6, alpha = 0.6) +
    geom_hline(yintercept = y_line, linetype = "dotted", color = "grey50") +
    simple_gradient +
    labs(x = "Time", y = "Effective N") +
    coord_cartesian(ylim = c(0, y_limit)) +
    base_theme
  if(time_scale == "log")(
    p1 <- p1 + scale_x_log10()
  )
  
  # --- Panel 2: Events ---
  p2 <- ggplot(df_summary, aes(x = events, y = n.eff, group = profile_id, color = surv)) +
    geom_line(linewidth = 0.6, alpha = 0.6) +
    geom_hline(yintercept = y_line, linetype = "dotted", color = "grey50") +
    simple_gradient +
    labs(x = "Number of Events", y = NULL) +
    coord_cartesian(ylim = c(0, y_limit)) +
    base_theme +
    theme(axis.text.y = element_blank())
  
  # --- Panel 3: Event Rate ---
  p3 <- ggplot(df_summary, aes(x = events_rate, y = n.eff, group = profile_id, color = surv)) +
    geom_line(linewidth = 0.6, alpha = 0.6) +
    geom_hline(yintercept = y_line, linetype = "dotted", color = "grey50") +
    simple_gradient +
    labs(x = "Proportion of Events", y = NULL) +
    coord_cartesian(ylim = c(0, y_limit)) +
    base_theme +
    theme(axis.text.y = element_blank())
  
  # --- Panel 4: Survival (The Convergence frontier) ---
  p4 <- ggplot(df_summary, aes(x = 1 - surv, y = n.eff, group = profile_id, color = surv)) +
    geom_line(linewidth = 0.4, alpha = 0.3) +
    # Highlight the max points frontier with a thicker black line
    geom_line(data = max_points, aes(x = 1 - surv, y = n.eff, group = profile_id), 
              color = "black", linewidth = 0.8, alpha = 0.8, inherit.aes = FALSE) +
    geom_hline(yintercept = y_line, linetype = "dotted", color = "grey50") +
    simple_gradient +
    labs(x = "1 - Survival (Incidence)", y = NULL) +
    coord_cartesian(ylim = c(0, y_limit)) +
    base_theme +
    theme(axis.text.y = element_blank())
  
  # Combine using patchwork
  combined <- (p1 | p2 | p3 | p4) + 
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Simplified Effective Sample Size Dynamics",
      subtitle = "All unique covariate profiles colored by survival probability"
    )
  
  return(combined)
}