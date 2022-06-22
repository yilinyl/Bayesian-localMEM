# Simulation for 2-stage local MEM
library("foreach")
library("partitions")
library("svMisc")
source("2_stage_anal.R")
source("local_mem_utils.R")
source("basket.R")

library(clinfun)

ph2simon(0.15, 0.45, 0.1/4, 0.2)  # alpha=0.1 / 4
# r1 n1 r  n EN(p0) PET(p0)
# Optimal  1  6 7 25  10.25  0.7765
# Minimax  1 10 5 16  12.73  0.5443

# Simon's exact probability
# r1 <- 1
# n1 <- 10
# n2 <- 6
# r <- 5
# 
# # Exact probability of rejecting H0 for Simon's two-stage design
# p_rej0 <- sapply(c(0.15, 0.45), function(p) {
#   1-(pbinom(r1, n1, p) + sum(sapply(c((r1+1):min(n1, r)), function(x) 
#   {dbinom(x, n1, p) * pbinom(r-x, n2, p)})))
# })
# 
# p_rej0

pmat <- rbind(c(0.15, 0.15, 0.15, 0.15),
              c(0.15, 0.15, 0.15, 0.45),
              c(0.15, 0.15, 0.45, 0.45),
              c(0.15, 0.45, 0.45, 0.45),
              c(0.45, 0.45, 0.45, 0.45))

B <- 4
n_priors <- 3
p0 <- pmat[1, ]
nsim <- 5000

part_prior_mat <- calc_priors(B, n_priors)
part <- part_prior_mat[, 1:B]
prior_mat <- part_prior_mat[, -(1:B)]

n1_vec <- rep(10, 4)
n2_vec <- rep(16, 4)
mem_2stg_sim_B4 <- local_mem_2stg(n1_vec, n2_vec, pmat, prior_mat, part, nsim=nsim)

save.image(paste("B", B, "_2stg_", paste(p0, collapse = "_"),
                 ".RData", sep = ""))
# dim(mem_2stg_sim_B4)

fwer_ctrl <- fwer_or_power_2tg(mem_2stg_sim_B4[1, 1, , ,], epsilon = 0.001, 
                            gam1_range = c(0.7, 0.999), gam2_range = c(0.7, 0.999), 
                            h0_idx = pmat[1, ] <= p0, tune=T)
# qf1 <- 0.985
# qf2 <- 0.900
tune_h1 <- fwer_or_power_2tg(mem_2stg_sim_B4[1, 5, , ,], gam_grid = fwer_ctrl[, 1:2],
                  h0_idx = pmat[5, ]<=p0, tune=F)

save.image(paste("B", B, "_2stg_tune", paste(p0, collapse = "_"),
                 ".RData", sep = ""))
# power_tune_2stg(mem_2stg_sim_B4[1, 5, , ,], gam_grid = matrix(fwer_ctrl[1, ], nrow=1), power_target = 0)
# fwer1(post_prob[, 1, ])
# dim(mem_2stg_sim_B4[1, 1, , ,])


