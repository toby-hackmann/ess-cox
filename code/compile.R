##### THIS SCRIPT COMPILES ALL REQUIRED FUNCTIONS FOR THE EFFECTIVE N PROJ #####

### Packages required ###
pkgs <- c(
  "survival",
  "tinytex",
  "tidyverse",
  "dplyr",
  "ggplot2",
  "plotly",
  "simsurv",
  "matrixStats",
  "ggsurvfit",
  "svglite",
  "zoo",
  "mstate",
  "patchwork",
  "grid",
  "knitr",
  "surveillance",
  "here"
)

options(repos = list(CRAN="http://cran.rstudio.com/"))

# Check if the package is installed, if yes load, if no: install + load 
vapply(pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  require(pkg, character.only = TRUE, quietly = TRUE)
}, FUN.VALUE = logical(length = 1L))

rm(pkgs)

colon_df <- colon
names <- c('sex', 'obstruct', 'perfor', 'differ', 'extent', 'surg', 'node4', 'etype')
colon_df <- colon_df |> 
  mutate( across( all_of(names), as.numeric ) )
rm(names)
colon_df <- colon_df[colon_df$etype == 2, ]

col <- c("lp" = "#9f84af", "dp" = "#37293F", "r" = "#C2666B", "y" = "#c6aa2c",  "b" = "#2E7691" )

source(here::here("code", "modified.R"))
source(here::here("code", "calculate_ess.R"))
source(here::here("code", "survfit_n.R"))
source(here::here("code", "plot_effective_n.R"))
source(here::here("code", "variance_cumhaz.R"))
source(here::here("code", "bootstrap.R"))
source(here::here("code", "variance_mstate.R"))
source(here::here("code", "sf_to_df.R"))
