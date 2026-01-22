
####################### EFFECTIVE N FOR SURVIVAL AT TIME #####################

calculate_ess <- function( df, 
                           n,
                           round = TRUE, 
                           uncensored = TRUE, 
                           bounds = TRUE,
                           mod = TRUE
){
  
  # Calculate the variance from the standard error and survival
  df$var <- df$surv^2*df$std.err^2
  
  ##### Calculate standard effective N #####
  n.eff <- (df$surv*(1-df$surv))/df$var
  
  # We set the effective sample size equal to the risk set if no event occurred
  n.eff[which( cumsum(df$n.event) == 0 )] <- df$n.risk[which( cumsum(df$n.event) == 0 )]
  
  # effective sample size is undefined when risk set or surv is 0
  # we set it to the last known value at that point
  if( any(is.nan(n.eff)) ) 
    n.eff[ which( is.nan(n.eff) ) ] <- n.eff[ min( which( is.nan(n.eff) ) ) - 1 ]
  
  # Rounding
  if( round ) n.eff <- round( n.eff )
  
  # Add to the output obj
  out <- list( n.eff = n.eff )
  
  
  ##### Calculate the modified effective sample size #####
  if( mod ){
    n.eff.mod <- modified( df )
    
    # We set the effective sample size equal to the risk set if no event occured
    n.eff.mod[which( cumsum(df$n.event) == 0 )] <- df$n.risk[which( cumsum(df$n.event) == 0 )]
    
    # Effective sample size is undefined when risk set or surv is 0
    # we set it to the last known value at that point
    if( any(is.nan(n.eff.mod)) )
      n.eff.mod[ which( is.nan(n.eff.mod) ) ] <- n.eff.mod[ min( which( is.nan(n.eff.mod) ) ) - 1 ]
    
    # Round
    if( round ) n.eff.mod <- round( n.eff.mod )
    out <- append(out, list( n.eff.mod = n.eff.mod ))
  }
  
  # Calculate number not censored
  if( uncensored ){
    n.uncensor <- n - cumsum( df$n.censor )
    out <- append( out, list( n.uncensor = n.uncensor))
  }
  
  #calculate lower and upper bounds
  if( bounds ){
    n.lower <- (df$n.risk-df$n.event)/df$surv
    n.upper <- cumsum(df$n.event)/(1-df$surv)
    out <- append( out, list( n.lower = n.lower, n.upper = n.upper))
  }
  return(out)
}
