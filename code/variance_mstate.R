##### CALCULATE THE VARIANCE USING MSTATE #####

variance_mstate <- function( cox, newdata ){
  
  # First, we define the transition matrix of a Cox model
  tmat <- transMat(list(c(2), c()), names = c("Alive", "Dead"))
  
  # Next, we prepare the data from Cox into the MSTATE format
  data <- data.frame( time = cox[["y"]][,1], status = cox[["y"]][,2], model.matrix(cox) )
  data.ms <<- msprep( time = c(NA, "time"), 
                  status = c( NA, "status"),
                  data = data, trans = tmat, 
                  keep = colnames(model.matrix(cox) ) )
  newdata$strata <- 1
  
  # Now, calculate the new cox object needed for msfit
  ms_cox <- coxph( update(cox$formula, reformulate(c( ".", "strata(trans)") ) ), data = data.ms )
  
  # Fit the multistate model
  ms_fit <- msfit( ms_cox, trans = tmat, newdata = newdata, vartype = "aalen" )
  rm(data.ms, envir = as.environment(.GlobalEnv) )
  
  # Extend ms_fit to the full time vector
  varhaz <- data.frame(time = ms_fit$varHaz$time, var = ms_fit$varHaz$varHaz )
  varhaz <- merge(varhaz, data.frame(time = sort(data$time)), by = c("time", "time"), all.y = TRUE)
  varhaz$var <- zoo::na.locf(varhaz$var, fromLast = FALSE)
  varhaz <- distinct(varhaz, time, .keep_all = TRUE)
  
  # Return the standard errors
  return( sqrt( varhaz$var ) )
}
