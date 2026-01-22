plot_effective_n <- function( 
    sf, 
    title = NA, 
    legend = c(0.15, 0.4),
    no_legend = FALSE,
    bounds = TRUE,
    survival = FALSE,
    max = FALSE,
    ylim = NA,
    mod = FALSE,
    xlim = NA,
    labels = TRUE,
    ticks_y = TRUE,
    ticks_x = TRUE,
    ticks_y_sec = TRUE,
    xmax = NA,
    xlab = "Time (days)",
    ylab = "Sample Size"
){
  
  if( class(sf) == "data.frame" ){
    df <- sf
    sf <- list( time = df$time,
                surv = df$surv,
                n.risk = df$n.risk,
                n.eff = df$n.eff,
                n = max(df$n.risk),
                n.uncensor = df$n.uncensor)
    if( mod ) sf <- append( sf, list( n.eff.mod = df$n.eff.mod))
    if( bounds ) sf <- append( sf, list( n.upper = df$n.upper, n.lower = df$n.lower))
    class(sf) <- "survfit"
  }
  
  
  if( length(sf$n) > 1 ){
    error( "Cannot plot a stratified Kaplan-Meier, use `plot_km_eff()` instead.")
  }
  
  # Define color values
  color_values <- c("Effective"= "#37293F", "At Risk" = "#C2666B", "Uncensored" = "#c6aa2c")
  if (bounds) {
    color_values <- c(color_values, "Bounds" = "#c6b5cf")
  }
  if (survival) {
    color_values <- c(color_values, "Survival" = "#2E7691")
  }
  if (mod) {
    color_values <- c(color_values, "Modified" = "#9f84af")
  }
  
  # Define line types
  linetype_values <- c("Effective"= "solid", "At Risk" = "dashed", "Uncensored" = "dotdash")
  if (bounds) {
    linetype_values <- c(linetype_values, "Bounds" = "dashed")
  }
  if (survival) {
    linetype_values <- c(linetype_values, "Survival" = "solid")
  }
  if (mod) {
    linetype_values <- c(linetype_values, "Modified" = "solid")
  }
  
  # X-axis labeling
  if( is.na(xmax) ){
    xmax <- ceiling(max(sf$time))
  }
  
  # Create the plot
  plot <- ggplot()
  if (bounds) {
    plot <- plot +
      geom_step(aes(x = sf$time, y = sf$n.upper, color = "Bounds", linetype = "Bounds"), linewidth = 1.3, alpha = 0.5) +
      geom_step(aes(x = sf$time, y = sf$n.lower, color = "Bounds", linetype = "Bounds"), linewidth = 1.3, alpha = 0.5)
  }
  
  if (survival){
    plot <- plot +
      geom_step(aes(x = sf$time, y = (1-sf$surv) * sf$n, color = "Survival", linetype = "Survival"), linewidth = 1.3) +
      scale_y_continuous(
        breaks = seq(0, min(2 * sf$n, max(c(sf$n.eff, sf$n))), length.out = 5),
        limits = c(0, min(2 * sf$n, max(c(sf$n.eff, sf$n))))
      ) +
      scale_y_continuous(sec.axis = sec_axis(~./sf$n, name = "Event probability"))
  }
  
  # Adds modified effective sample size
  if (mod) {
    plot <- plot +
      geom_step(aes(x = sf$time, y = sf$n.eff.mod, color = "Modified", linetype = "Modified"), linewidth = 1.3)
  }
  
  # Adds maximum line
  if ( max ) {
    plot <- plot + 
      geom_hline(yintercept = sf$n, linetype = "dashed", color = "black", linewidth = 1.3)
  }
  
  plot <- plot +
    geom_step(aes(x = sf$time, y = sf$n.eff, color = "Effective", linetype = "Effective"), linewidth = 1.3) +
    geom_step(aes(x = sf$time, y = sf$n.risk, color = "At Risk", linetype = "At Risk"), linewidth = 1.3) +
    geom_step(aes(x = sf$time, y = sf$n.uncensor, color = "Uncensored", linetype = "Uncensored"), linewidth = 1.3) +
    labs(x = xlab, y = ylab) +
    theme_minimal() +
    theme(legend.position.inside = legend, panel.grid.major.x = element_blank()) +
    theme(text = element_text(size = 12), plot.title = element_text( size = 18, hjust = 0.3)) +
    scale_color_manual(values = color_values, name = "Legend") +
    scale_linetype_manual(values = linetype_values, name = "Legend")
  
  if( length(unique(sf$time)) > 5 ){
    plot <- plot + scale_x_continuous(breaks = seq(0, xmax, length.out = 5))
  }
  
  # Only include title if it is given
  if( !is.na(title) ){
    plot <- plot +
      labs( title = title )
  }
  
  # Removes legend
  if( no_legend ){
    plot <- plot +
      theme(legend.position = "none")
  }
  
  # Removes axis labels
  if( !labels ){
    plot <- plot +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
            plot.tile = element_blank())
  }
  
  # Overrides default ylims
  if( all(!is.na(ylim)) ){
    plot <- plot +
      scale_y_continuous( limits = ylim )
  }
  # Overrides default xlims
  if( all(!is.na(xlim)) ){
    plot <- plot +
      scale_x_continuous( limits = xlim )
  }
  
  # Removes axis ticks and numbers
  if( !ticks_y ){
    plot <- plot +
      theme( axis.text.y.left = element_blank() )
  }
  if( !ticks_x ){
    plot <- plot +
      theme( axis.text.x = element_blank() )
  }
  if( !ticks_y_sec ){
    plot <- plot +
      theme( axis.text.y.right = element_blank() )
  }
  
  return(plot)
}
