extract_sf <- function(sf_obj, tx, prevalence){
  df <- sf_to_df(sf_obj, time_grid)
  df$tx <- factor(tx)
  df$prevalence <- prevalence
  df
}

run_one <- function(rep_id, scale, shape, maxtime, hr, n, time_grid,
                    model_type = c("standard", "robust")){
  
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
