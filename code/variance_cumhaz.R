
####################### COX STD ERR - KLEIN MOESCHBERGER #####################
# This script starts with some simple supporting functions that are used in the
# Klein & Moeschberger variance calculations in such a way that it resembles
# the method from the book closely


############# SUPPORTING FUNCTIONS ####################

# Calculate the weighted risk set
.W <- function( t, cox, derivative = FALSE ){
  df <- data.frame( time = cox[["y"]][,1], lin.pred = cox[["linear.predictors"]] )
  
  if( derivative ){
    X <- model.matrix( cox )
    Y <- numeric( ncol(X) )
    
    for( j in 1:ncol(X) ){
      Y[j] <- ( X[ df[, 1] >= t, j]*exp( df[ df[, 1] >= t, "lin.pred" ] - cox$mean.risk ) ) %>% sum()
    }
    return( Y )
  }
  else sum( exp( df[ df[, 1] >= t, "lin.pred" ] - cox$mean.risk ) )
}


# Risk of actual events <- shitty idea
.R <- function( t, cox ){
  # 0 Risk if there are no events
  if( sum( cox$y[ , 2][cox$y[, 1] ==  t ] ) == 0 ) return( 0 )
  
  # Otherwise, calculate the total risk of the events at that time
  return( sum( exp( cox$linear.predictors[ cox$y[, 1] == t & cox$y[, 2] == 1] ) ) )
}


# Variance of baseline cumulative hazard, with choice
.Q1 <- function( df, method, new_risk ){

  # For Klein, we calculate the basic Aalen variance
  if( method %in% c("Klein") ) return( cumsum( df$n.event/(df$risk^2 ) ) ) 
  
  # For Tsiatis, we calculate the Aalen variance multiplied by squared HR
  else if( method %in% c("Tsiatis") ) return( new_risk^2*cumsum( df$n.event/(df$risk)^2 ) )
      
  # If the method is Link, we calculate the p_i, and from this the GW product    
  else if( method %in% c("Link") ){
    p_i <- exp( - (new_risk*df$n.event)/df$risk )
    return( (cumprod( (1 + (1-p_i)/(df$n.risk*p_i) ) ) - 1 ) )
    #return( cumsum( (1-p_i)/(df$n.risk*p_i) ) )
  }
  
  else stop( "Baseline not recognized")
}


# Supporting function for Q2
.Q3 <- function( df, cox, newdata, method, new_risk ){
  
  # First, we need the derivative vector over time in all directions
  M <- lapply( df$time, .W, cox, TRUE ) %>% do.call( rbind, . )
  
  # Now we loop over the variables in the newdata and calculate the vector
  # element each time using the Klein & Moeschberger formula
  if( method %in% c("Klein") ){
    q3 <- matrix(0, nrow = nrow(M), ncol = ncol(newdata))
    for( k in 1:ncol(newdata)){
      q3[ , k] <- cumsum( (M[ , k]/df$risk - newdata[ , k])*(df$n.event/df$risk) )
    }
  }
  
  # Formula for Link, with additional new_risk terms
  else if ( method %in% c( "Link" ) ){
    q3 <- matrix(0, nrow = nrow(M), ncol = ncol(newdata))
    for( k in 1:ncol(newdata)){
      # I grayed out the newdata term - this is not supposed to be here, but it does make the slope of the two groups the same - something to be said for that? 
      q3[ , k] <- new_risk * cumsum( (M[ , k]*df$n.event/df$risk^2) ) - cumsum(df$n.event/df$risk) * new_risk * newdata[ , k] 
    }
  }
  
  # Formula for Tsiatis, different form, but should be same as Klein
  else if ( method %in% c( "Tsiatis" ) ){
    q3 <- matrix(0, nrow = nrow(M), ncol = ncol(newdata))
    for( k in 1:ncol(newdata)){
      q3[ , k] <- cumsum( (M[ , k]*df$n.event/df$risk^2) ) - cumsum(df$n.event/df$risk) * newdata[ , k]
    }
  }
  
  else stop( "Method not recognized")
  
  return(q3)
}


# Variance depending on the coefficient estimation
.Q2 <- function( df, cox, newdata, method, new_risk ){
  
  # Since Klein & Moeschberger don't specify how to calculate the variance of
  # the coefficient estimates we just take the standard errors from the cox
  # model as calculated by the survival package.
  
  # Tsiatis, according to Link, multiplies by squared risk in this step
  if( method == "Tsiatis" ) return( new_risk^2 * rowSums( (.Q3( df, cox, newdata, method, new_risk ) %*% cox$var) * .Q3( df, cox, newdata, method, new_risk ) ) )
  
  # Otherwise
  rowSums( (.Q3( df, cox, newdata, method, new_risk ) %*% cox$var) * .Q3( df, cox, newdata, method, new_risk ) )
}



############ ACTUAL KLEIN & MOESCHBERGER FUNCTION ##################
variance_cumhaz <- function( df, cox, newdata, method = "Tsiatis", which = "default" ){
  # Where df is the dataframe of the survfit object, cox the cox model and
  # X the design matrix <---- this last part needs to be resolved before this
  # method works on new data
  
  # First we add the value of the weighted risk set to the dataframe
  #cox$mean.risk <- cox$means * cox$coefficients %>% sum()
  cox$mean.risk <- mean(cox$linear.predictors)
  df$risk <- lapply( df$time, .W, cox ) %>% unlist()
  new_risk <- exp( sum( newdata * cox$coefficients ) - cox$mean.risk )
  
  Q1 <- .Q1( df, method, new_risk )
  Q2 <- .Q2( df, cox, newdata, method, new_risk )
  Q3 <- Q1 + Q2
  
  return(list(std.err = sqrt(Q3), 
              std.chaz = sqrt(Q1), 
              std.coef = sqrt(Q2))
  )
}


