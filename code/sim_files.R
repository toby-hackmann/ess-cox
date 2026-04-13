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