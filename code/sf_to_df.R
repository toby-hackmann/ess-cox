
##### TURN A SURVFIT OBJECT INTO A DATAFRAME WITH A CERTAIN TIME SCALE #####

sf_to_df <- function( sf, time ){
  # Turn list in to df
  df <- unclass(sf)[c("time", "surv", "lower", "upper", "n.risk", "n.uncensor", "n.eff", "n.lower", "n.upper", "std.err", "std.chaz")] |>  
    data.frame()
  
  if( "std.coef" %in% names(sf) ){
    df.temp <- unclass(sf)["std.coef"] |> data.frame()
    df$std.coef <- df.temp$std.coef
  }
  
  # merge with time grid
  df <- merge( df, data.frame(time = time ), by = "time", all = TRUE, sort = TRUE )
  
  # Fill the dataframe downwards
  df <- fill(df, 2:8, .direction = "downup")
  
  return(df[ df$time %in% time, ])
  
}