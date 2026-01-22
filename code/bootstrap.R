########################## COX STD ERR - BOOTSTRAP #############################
# This script uses Bootstrapping to calculate the standard errors of the 
# cumulative hazard from a Cox model. These will be used in the calculation of
# effective sample size for a prediction using the Cox model.


### BOOTSTRAPPING PROBABLY DOES NOT WORK

# First, we write the short function of what is done every bootstrap iteration
iteration <- function( df, cox, newdata, time ){
  
  # Step 1: Fit the cox model on the bootstrapped dataframe
  model <- coxph( cox$formula, data = df )
  
  # Step 2: Fit on the newdata
  sf <- survfit( model, newdata = newdata )
  
  # Step 3: Get the cumulative hazard
  cumhaz <- summary(sf, times = time)$cumhaz
  
  # Step 4: Enlargen the cumulative hazard until it is length time by adding 0s
  cumhaz <- c( rep(0, length(time)-length(cumhaz)), cumhaz )
  
  # Step 5: Return the cumulative hazard
  return( cumhaz )
}

bootstrap <- function( cox, newdata, B ){
  
  # First, we need the data from the cox model
  df <- data.frame( time = cox[["y"]][,1], status = cox[["y"]][,2], model.matrix(cox) )
  
  # We need to get the vector of times where we want to evaluate the cumhaz
  time <- unique(df$time)
  
  # Now, for each bootstrap iteration, we need to return the predicted cumhaz
  # at all times in `time`
  cumhaz <- replicate( B, iteration( df[sample(nrow(df), replace = TRUE),], cox, newdata, time ) )
  
  # Now calculate the standard error by row
  se <- apply( cumhaz, 1, sd )
  
  # Return the standard error
  return(se)
}


