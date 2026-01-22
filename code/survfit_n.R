
####################### EFFECTIVE N FOR SURVIVAL AT TIME #####################

survfit_n <- function( sf, 
                       cox = NULL,
                       round = TRUE, 
                       uncensored = TRUE, 
                       bounds = TRUE,
                       mod = TRUE,
                       method = "Default",
                       B = 100, 
                       time = NULL,
                       coef = F,
                       chaz = F
){
  # Check if, when the sf object is a survfit prediction of a Cox model, the cox
  # model is also provided to the function
  if( all("survfitcox" %in% class(sf), class(cox) != "coxph") ){
    error( "Providing the coxph object is required when looking at effective
           sample size for a Cox model prediction.")
  }
  
  # Validity checks of method
  if( all(!method %in% c("Kaplan-Meier"), is.null(cox), !"survfitcox" %in% class(sf) ) ){
    warning( "The survfit object is not a Cox model, so the method is defaulted back to Kaplan-Meier" )
    method <- "Kaplan-Meier"
  }
  
  if( all( !method %in% c("Link", "Default", "Klein", "Bootstrap", "Tsiatis", "MSTATE"), "survfitcox" %in% class(sf) ) ){
    error( "Method must be one of 'Default', 'Link', 'Klein', 'Tsiatis', 'MSTATE' or 'Bootstrap'." )
  }
  
  # First, we unwrap the survfit object into a dataframe with all relevant data
  if( any( nrow(sf$newdata) == 1, !("survfitcox" %in% class(sf))) )  
    df <- unclass(sf)[c("time", "surv", "n.risk", "n.event", "n.censor", "std.err", "std.chaz")] |>  
      data.frame()
  else if( !is.null(time) ){
    index <- which.min(sf$time < time)
    df <- data.frame(surv = unclass(sf$surv[index, ]), std.err = unclass(sf$std.err[index, ]))
    return( (1-df$surv)/(df$surv*df$std.err^2) )
  }
  else
    error("You can only get the effective sample size for one new patient, or have to specify at what time to evaluate effective sample size for many patients.")
  
  # If it is a cox fit, we need to calculate the stabilizer - not sure if this
  # makes sense, just comment out for now. This scales it to the complete sample
  #if( "survfitcox" %in% class(sf) ){
  #  # Calculate the weighted risk set at each time
  #  w.risk <- NULL
  #  for( i in 1:length(df$time) ){
  #    index <- which( unclass(cox$y)[, 1] >= df$time[i] ) 
  #    w.risk[i] <- sum( exp( cox$linear.predictors[index] ) )
  #  }
  #  stabiliser <- df$n.risk / w.risk
  #  # Multiply by individual linear predictor
  #  stabiliser <- stabiliser * as.numeric(exp(cox$coefficients*cox_fit$newdata))
  #} else stabiliser = 1
  
  # If it is a cox fit, we may want to calculate the std errors differently from
  # the default survival package method
  if( "survfitcox" %in% class(sf) & (coef | chaz) ){
    mod <- FALSE
    newdata <- model.matrix(cox, sf$newdata)
    out <- switch( method,
      "Default" = variance_cumhaz( df, cox, newdata, "Tsiatis" ),
      "Link" = variance_cumhaz( df, cox, newdata, method ),
      "Klein" = variance_cumhaz( df, cox, newdata, method ),
      "Tsiatis" = variance_cumhaz( df, cox, newdata, method ),
      "Bootstrap" = bootstrap( cox, newdata, B),
      "MSTATE" = variance_mstate( cox, newdata )
    )
    df$std.err <- out$std.err
    sf$std.err <- df$std.err
    if( chaz ){
      sf$std.chaz <- out$std.chaz
    }
    if( coef ){
      sf$std.coef <- out$std.coef
    }
  }



  # Regular Kaplan-Meier without strata
  if( length(sf$n) == 1 & !("survfitcox" %in% class(sf)) ){
    ess <- calculate_ess( df, sf$n, 
                          round = round, 
                          uncensored = uncensored, 
                          bounds = bounds, 
                          mod = mod)
    #ess[["n.eff"]] <- ess[["n.eff"]]*stabiliser
    #ess[["n.eff.mod"]] <- ess[["n.eff.mod"]]*stabiliser
    #sf$std.chaz <- sqrt(cumsum(sf$n.event/sf$n.risk^2))
    sf <- append( sf, ess )
    class(sf) <- "survfit"
    return(sf)
  }
  # Multiple strata
  else if( !("survfitcox" %in% class(sf))){
    # Find the cut points of the different strata
    cut <- c(0, cumsum(sf$strata) )
    # Run once outside loop for list initializaion
    ess <- calculate_ess( df[(cut[1]+1):cut[2], ], sf$n[1], 
                          round = round, 
                          uncensored = uncensored, 
                          bounds = bounds, 
                          mod = mod )
    for( i in 2:length(sf$strata) ){
      ess <- Map(c, ess, calculate_ess( df[(cut[i]+1):cut[i+1], ], sf$n[i], 
                                        round = round, 
                                        uncensored = uncensored, 
                                        bounds = bounds, 
                                        mod = mod)) 
    }
    sf <- append( sf, ess )
    class(sf) <- "survfit"
    return( sf )      
  }
  else{
    ess <- calculate_ess( df, sf$n, 
                          round = round, 
                          uncensored = uncensored, 
                          bounds = bounds, 
                          mod = mod)
    sf <- append( sf, ess )
    class(sf) <- "survfit"
    return( sf )
  }
  class(sf) <- "survfit"
  return( sf )
}
