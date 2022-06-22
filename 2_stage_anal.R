# updated - 2 stage analysis

local_mem_2stg <- function(n1, n2, pmat, prior_mat, all_part,
                           B=ncol(pmat), n_prior=ncol(prior_mat), 
                           a0=1, b0=1, nsim=5000, seed=12321) {
  # mem_pp: (n_prior, n_scenario, n_sim, n_stage, n_basket)
  mem_pp <- array(NA, dim=c(n_prior, nrow(pmat), nsim, 2, B))
  p0 <- pmat[1, ]
  set.seed(seed)
  for(prior_choice in 1:n_prior){
    prior_part <- prior_mat[, prior_choice]
    #given a prior choice, run through scenarios
    for(j in 1:nrow(pmat)){
      ptrue <- pmat[j, ]
      x1 <- sapply(1:nsim, function(i) {rbinom(B, size = n1, prob = ptrue)})
      x2 <- x1 + sapply(1:nsim, function(i) {rbinom(B, size = n2-n1, prob = ptrue)})
      
      for(i in 1:nsim){
        # 2-stage analysis
        # stage 1
        mem_pp[prior_choice, j, i, 1, ] <-
          local_mem_unit(x1[, i], n1, ptrue, p0, B, prior_part, all_part,
                         a0=a0, b0=b0)
        # stage 2
        mem_pp[prior_choice, j, i, 2, ] <-
          local_mem_unit(x2[, i], n2, ptrue, p0, B, prior_part, all_part,
                         a0=a0, b0=b0)
        
      }
      if(j==nrow(pmat)) message("DONE!")
    }
  }
  return(mem_pp)
}


# local MEM for single stage (analysis unit)

local_mem_unit <- function(x, Nvec, ptrue, p0, B, prior_part, part, 
                           a0=1, b0=1) {
  # x <- rbinom(B, size = Nvec, prob = ptrue)
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
  # mem_pp[prior_choice, j, i, ] <- 1 - pbeta(p0, a1, b1)
  return(1 - pbeta(p0, a1, b1))
  # print(c(i, j, prior_choice))
}

