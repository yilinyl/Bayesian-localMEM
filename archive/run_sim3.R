# Example code for running simulation (for 2-stage analysis)

if (!require(partitions)) {
  install.packages("partitions")
  library(partitions)
}
if (!require(pscl)) {
  install.packages("pscl")
  library(pscl)
}

source("./sequential_analysis.R")
source("./basket_part.R")
source("./basket_cluster.R")
source("./summary_basket.R")
source("./utils.R")

pmat <- rbind(c(0.15, 0.15, 0.15, 0.15, 0.15, 0.15),
              c(0.15, 0.15, 0.15, 0.15, 0.15, 0.45),
              c(0.15, 0.45, 0.15, 0.15, 0.15, 0.45),
              c(0.45, 0.45, 0.15, 0.15, 0.15, 0.45),
              c(0.45, 0.45, 0.45, 0.15, 0.15, 0.45),
              c(0.45, 0.45, 0.45, 0.45, 0.15, 0.45),
              c(0.45, 0.45, 0.45, 0.45, 0.45, 0.45),
              c(0.45, 0.45, 0.15, 0.35, 0.35, 0.45))

nsims <- 5000
# int_n <- c(10, 9)
Nmax <- rep(19, 6)
Nmin = 10
n_basket <- length(Nmax)
props <- rep(1 / n_basket, n_basket)

sim3_new <- lapply(1:nrow(pmat), function(i) {
  out <- lapply(1:nsims, function(s) {
    if (s %% 100 == 0) {
      print(paste("---- Simulation", s, "/", nsims, "----"))
    }
    sequential_monitor(pmat[i, ], props, Nmin, Nmax_v = Nmax, 
                       delta = 1,  bmin = n_basket, step = 9, verbose = T,
                       stepwise = TRUE, evaluate = FALSE)
  })
})

# save(sim3_new, file = "./sim3_new1.RData")


