rm(list=ls(all=TRUE))
library("foreach")
library("partitions")
library("svMisc")

source("basket.R")
# ---------------------------------------------------------
# set up the basket trial
# ---------------------------------------------------------
# the number of baskets
B <- 4
# Nvec: the sample size of each basket
Nvec <- rep(19, B)

# enumerate all possible partitions of the B baskets
part <- get.part(R = B, max_cl = B)
K <- nrow(part)

# number of blocks/unique response rate in each partition
n_bk <- apply(part, 1, function(x){
  #max(table(x))-1
  length(unique(x))
})

# the response scenarios for simulation
# each row of pmat is a different simulation scenario

# pmat <- rbind(c(0.15, 0.15, 0.15, 0.15, 0.15, 0.15),
#               c(0.15, 0.15, 0.15, 0.15, 0.15, 0.45),
#               c(0.15, 0.15, 0.15, 0.15, 0.45, 0.45),
#               c(0.15, 0.15, 0.15, 0.45, 0.45, 0.45),
#               c(0.15, 0.15, 0.45, 0.45, 0.45, 0.45),
#               c(0.15, 0.45, 0.45, 0.45, 0.45, 0.45),
#               c(0.45, 0.45, 0.45, 0.45, 0.45, 0.45))

# pmat <- pmat[1:5, 3:6]

# pmat <- rbind(c(0.25, 0.25, 0.15, 0.15),
#               c(0.25, 0.25, 0.15, 0.45),
#               c(0.25, 0.25, 0.45, 0.45),
#               c(0.25, 0.55, 0.45, 0.45),
#              c(0.55, 0.55, 0.45, 0.45))

# fixed response rates for all baskets: 
#     0.15 under the null vs. 0.45 under the alternative
pmat <- rbind(c(0.15, 0.15, 0.15, 0.15),
              c(0.15, 0.15, 0.15, 0.45),
              c(0.15, 0.15, 0.45, 0.45),
              c(0.15, 0.45, 0.45, 0.45),
              c(0.45, 0.45, 0.45, 0.45))
print(pmat)

# null response rate
p0 <- pmat[1, ]
# beta prior for each basket
a0 <- b0 <- 1
# ---------------------------------------------------------
# set up simulation
# ---------------------------------------------------------
#number of simulations
nsim <- 5000

# array for saving the posterior prob:
# the 1st index for number of mem priors
# the 2ed index for number of scenarios
# the 3rd index for number of simulations
# the 4th index for number of baskets
mem_pp <- array(NA, dim=c(3, nrow(pmat), nsim, B))

prior_mat <- matrix(NA, nrow = K, ncol = 3)
# set the prior choice
for(prior_choice in 1:3){
  # has no preference on which partition to choose
  if(prior_choice==1){
    prior_mat[, prior_choice] <- rep(1/K, K)
  }
  # promote local borrowing, believes in heterogeneity
  if(prior_choice==2){
    prior_mat[, prior_choice] <- n_bk/sum(n_bk)
    #prior_mat[, prior_choice] <- n_bk/sum(n_bk)
  }
  # has no preference on the number of blocks
  if(prior_choice==3){
    prior_mat[, prior_choice] <- n_bk^2/sum(n_bk^2)
    # this prior gives too much weight on naieve pooling and no pooling
    # and too little weight to each of the remaining partitions
    # for(k in 1:K){
    #   prior_mat[k, prior_choice] <- (1/B)/sum(n_bk==n_bk[k])
    # }
  }
}

colSums(prior_mat)
aggregate(prior_mat[, 1], by = list(n_bk), sum)
aggregate(prior_mat[, 2], by = list(n_bk), sum)
aggregate(prior_mat[, 3], by = list(n_bk), sum)

for(prior_choice in 1:3){
  prior_part <- prior_mat[, prior_choice]
  #given a prior choice, run through scenarios
  for(j in 1:nrow(pmat)){
    ptrue <- pmat[j, ]
    
    for(i in 1:nsim){
      set.seed(i)
      x <- rbinom(B, size = Nvec, prob = ptrue)
      y <- Nvec - x
      
      post_part <- update.part(x = x, n = Nvec, prior_part = prior_part, 
                               part = part)
      # partition with the maximum pp
      part_hat <- unlist(part[which.max(post_part), ])
      p_hat <- post_part[which.max(post_part)]
      
      a1 <- b1 <- NULL
      for(b in 1:B){
        # exchangeability of B baskets
        exc <- rep(0, B)
        exc[part_hat==part_hat[b]] <- p_hat
        exc[b] <- 1
        # print(exc)
        a1[b] <- a0 + sum(x*exc)
        b1[b] <- b0 + sum(y*exc)
      }
      mem_pp[prior_choice, j, i, ] <- 1 - pbeta(p0, a1, b1)
      
      print(c(i, j, prior_choice))
    }
    if(j==nrow(pmat)) message("DONE!")
  }
}

save.image(paste("./data/B", B, "_", paste(p0, collapse = "_"),
                   ".RData", sep = ""))
