# ------ Simulation --------
# run basket simulations
# setwd("~/Documents/Yale/Research/clinical_trial/script/basket_trial_local")

if (!require(partitions)) {
  install.packages("partitions")
  library(partitions)
}
if (!require(pscl)) {
  install.packages("pscl")
  library(pscl)
}

# source("./sequential_analysis.R")
# source("./basket_part.R")
# source("./basket_cluster.R")
# source("./summary_basket.R")
# source("./utils.R")

# ---- monitoring a basket trial ----
#' @param pmat probability matrix for simulation scenarios (#scenarios x #baskets)
#' @param props proportion of each subtype [vector]
#' @param nsims number of simulations
#' @param Nmin minimum sample size for a basket to enter interim analysis
#' @param Nmax maximum sample size for each basket to finish interim analysis [vector]
#' @param bmin minimum number of basket required to enter first interim analysis
#' @param stepwise if TRUE, then all baskets enter interim analysis with the same interim size [default=TRUE]
#' @param step an interim look is conducted each time when additional `step` patients are enrolled (only used when `stepwise`=TRUE)
#' @param evaluate if FALSE, baskets will be considered until reach maximum sample size (no early termination) [default=FALSE]
run_simulation <- function(pmat, props, nsims, Nmin, Nmax, bmin, delta=0.02,
                           stepwise=TRUE, step=1, evaluate=FALSE, verbose=FALSE, ...) {
  if (!is.null(pmat) & is.null(dim(pmat))) {
    pmat <- matrix(pmat, nrow = 1)
  }
  sim_data <- lapply(1:nrow(pmat), function(i) {
    out <- lapply(1:nsims, function(s) {
      if (s %% min(100, nsims) == 0) {
        print(paste("---- Simulation", s, "/", nsims, "----"))
      }
      sequential_monitor(pmat[i, ], props, Nmin, Nmax_v = Nmax, 
                         delta = delta,  bmin = bmin, step = step, verbose = verbose,
                         stepwise = stepwise, evaluate = evaluate)
    })
  })
  return(sim_data)
}


# Example:
# pmat <- rbind(c(0.15, 0.15, 0.15, 0.15, 0.15, 0.15),
#               c(0.15, 0.15, 0.15, 0.15, 0.15, 0.45),
#               c(0.15, 0.45, 0.15, 0.15, 0.15, 0.45),
#               c(0.45, 0.45, 0.15, 0.15, 0.15, 0.45),
#               c(0.45, 0.45, 0.45, 0.15, 0.15, 0.45),
#               c(0.45, 0.45, 0.45, 0.45, 0.15, 0.45),
#               c(0.45, 0.45, 0.45, 0.45, 0.45, 0.45),
#               c(0.45, 0.45, 0.15, 0.35, 0.35, 0.45))
# #Case 1: Balanced sample size, one-stage analysis
# Nmax <- rep(19, 6) 
# Nmin = Nmax
# n_basket <- length(Nmax)
# props <- rep(1 / n_basket, n_basket)
# nsims <- 5
# 
# run_simulation(pmat[1, ], props, nsims, Nmin, Nmax, bmin = n_basket, verbose = T)

# #Case 3: Balanced sample size, two-stage analysis n=(10, 19)
# Nmax <- rep(19, 6) 
# Nmin = 10
# n_basket <- length(Nmax)
# props <- rep(1 / n_basket, n_basket)
# nsims <- 5
# 
# run_simulation(pmat[1, ], props, nsims, Nmin, Nmax, bmin = n_basket, 
#                stepwise=T, step=9, verbose = T)



