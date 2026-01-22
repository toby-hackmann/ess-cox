
##################### MODIFIED EFFECTIVE SAMPLE SIZE ###########################
# This function slots in to the function survfit_n and is used to calculate the

modified <- function( df ){
  
  # Start with same values as normal
  surv.mod <- df$surv
  var.mod <- df$var
  
  # Change the survival and variances for those times when there are only censor
  for( i in which( df$n.event == 0 ) ){
    
    # Don't calculate before first event
    if( sum( df$n.event[1:i] ) == 0 ) next
    
    # When was the last event before the current censoring-only time?
    last <- max( which( df$n.event[1:i] != 0 ) )
    if( last == 1 ) next
    
    # Recalculate the survival for the current time
    surv.mod[i] <- df$surv[last - 1]*
                   (1 - ( df$n.event[last] - 1 )/df$n.risk[last] )*
                   (1 - 1/( df$n.risk[i] + 1 ) )
    
    # Recalculate the variance based for the current time
    var.mod[i] <- surv.mod[i]^2*( 
                    df$std.err[last - 1]^2 +
                    (df$n.event[last] - 1) /
                    (df$n.risk[last]*(df$n.risk[last] - df$n.event[last] + 1)) +
                    1/((df$n.risk[i]+1)*df$n.risk[i])
                    )
  }
  
  # Return the modified effective sample size
  return( surv.mod*(1-surv.mod)/var.mod )
}