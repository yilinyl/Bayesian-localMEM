#enumerate all possible partitions of baskets in a basket trial 
#and calculate the uncertainties related to each partition

##x: vector of responses in each basket
##n: vector of sample sizes in each basket
##beta prior for response rate: a0, b0 [default: a0 = b0 = 1]
##max_cl: maximum number of clusters a basket trial should be partitioned to

require("partitions")

# R: the number of trials to consider
# max_cl: maximum number of blocks
# same results as previous version
get.part <- function(R, max_cl){
  ## generate all the possible partitions 
  ## and store them in a matrix 
  part_mat <- t(setparts(R))
  part_mat <- part_mat[apply(part_mat, 1, function(x){
    length(unique(x))<=max_cl
  }), ]
  part_mat <- data.frame(part_mat)
  names(part_mat) <- LETTERS[1:R]
  #part_mat[-1, ]
  part_mat
}

update.part <- function(x, n, prior_part, part, a0 = 1, b0 = 1){
  R <- length(x)
  K <- nrow(part)
  
  p <- foreach(k = 1:K,.combine = "c")%do%{
    grp <- unlist(part[k, ])
    S <- aggregate(x, by=list(grp), sum)$x
    N <- aggregate(n, by=list(grp), sum)$x
    #calculate marginal probs m(s_j)
    prod((beta(a0+S, b0+N-S)/beta(a0,b0)))*prior_part[k]
  }
  #marginal probs m(s_j)
  mp <- sum(p)
  ## calculate posterior prob of each grouping structure
  p/mp
}

# post_prob: an nsim*B matrix of posterior probabilities
# nsim: number of simulated trials
# B: number of baskets
# gam: declare success in a basket if post_prob>gam
fwer1 <- function(post_prob, fwer_target = 0.10, epsilon = 0.001, gam = NULL){
  if(is.null(gam)){
    FWER <- 0
    gam <- 1
    while(FWER<fwer_target){
      gam <- gam - epsilon
      FWER <- mean(apply(post_prob>gam, 1, any))
    }
    gam <- gam + epsilon
    FWER <- mean(apply(post_prob>gam, 1, any))
  }else{
    FWER <- mean(apply(post_prob>gam, 1, any))
  }
  return(list("gam"=gam, "FWER"=FWER))
}

# tune posterior cutoff that controls FWER for 2-stage design
fwer_ctrl_2stg <- function(post_prob, fwer_target=0.1, epsilon = 0.001, 
                     gam1_range = c(0.7, 0.9), gam2_range=c(0.95, 1), h0_idx=NULL) {
  gam_grid <- expand.grid(seq(gam1_range[1], gam1_range[2], epsilon), 
                          seq(gam2_range[1], gam2_range[2], epsilon))
  print(paste("----- Test ", nrow(gam_grid), "combinations -----"))
  n_basket <- dim(post_prob)[3]
  if (is.null(h0_idx)) {
    h0_idx <- rep(TRUE, n_basket)
  }
  selected <- NULL
  for (i in 1:nrow(gam_grid)) {
    qf1 <- gam_grid[i, 1]
    qf2 <- gam_grid[i, 2]
    is_fut1 <- post_prob[, 1, ] <= qf1
    is_fut2 <- post_prob[, 2, ] <= qf2
    # is_eff <- !is_fut1 & !is_fut2
    is_fut <- is_fut1 | is_fut2
    fwer <- mean(sapply(1:nrow(is_fut), function(s) {any(!is_fut[s, h0_idx])}))
    print(fwer)
    if (fwer < fwer_target) {
      selected <- rbind(selected, c(qf1, qf2, fwer))
    }
  }
  return(selected)
}

# tune posterior cutoff that maximize power for 2-stage design
fwer_or_power_2tg <- function(post_prob, fwer_target=0.1, epsilon = 0.001, 
                              gam1_range = c(0.7, 0.9), gam2_range=c(0.95, 1),
                              gam_grid=NULL, h0_idx=NULL, tune=FALSE) {
  if (is.null(gam_grid)) {
    gam_grid <- expand.grid(seq(gam1_range[1], gam1_range[2], epsilon), 
                            seq(gam2_range[1], gam2_range[2], epsilon))
  }
  # gam_grid <- expand.grid(seq(gam1_range[1], gam1_range[2], epsilon), 
  #                         seq(gam2_range[1], gam2_range[2], epsilon))
  print(paste("----- Test ", nrow(gam_grid), "combinations -----"))
  n_basket <- dim(post_prob)[3]
  if (is.null(h0_idx)) {
    h0_idx <- rep(TRUE, n_basket)
  }
  h1_idx <- !h0_idx
  selected <- NULL
  for (i in 1:nrow(gam_grid)) {
    qf1 <- gam_grid[i, 1]
    qf2 <- gam_grid[i, 2]
    is_eff1 <- post_prob[, 1, ] > qf1
    is_eff2 <- post_prob[, 2, ] > qf2
    is_eff <- is_eff1 & is_eff2
    fwer <- mean(sapply(1:nrow(is_eff), function(s) {any(is_eff[s, h0_idx])}))
    err_power <- colMeans(is_eff)  # basket_wise_power or type I error
    power <- mean(err_power)
    if (tune) {  # tuning posterior cutoffs
      if (fwer < fwer_target) {
        selected <- rbind(selected, c(qf1, qf2, err_power, fwer, power))
      }
    }
    else {
      selected <- rbind(selected, c(qf1, qf2, err_power, fwer, power))
    }
  }
  if (!is.null(selected)) {
    colnames(selected) <- c("qf1", "qf2", LETTERS[1:n_basket], "FWER", "power")
  }
  return(selected)
}

# power for local mem
mem_eval <- function(mem_pp, pmat, B=NULL, n_prior=NULL, nsim=NULL, 
                     target=0.1) {
  if (is.null(B)) {
    B <- dim(mem_pp)[4]
  }
  if (is.null(n_prior)) {
    n_prior <- dim(mem_pp)[1]
  }
  if (is.null(nsim)) {
    nsim <- dim(mem_pp)[3]
  }
  p0 <- pmat[1, ]
  mem_pow <- array(NA, dim = c(n_prior, nrow(pmat), B))
  gam <- FWER <- NULL
  for(prior_choice in c(1:n_prior)){
    gam[prior_choice] <- fwer1(post_prob = mem_pp[prior_choice, 1, , ], 
                               fwer_target = target)$gam
    FWER[prior_choice] <- fwer1(post_prob = mem_pp[prior_choice, 1, , ], 
                                fwer_target = target)$FWER
  }
  print(paste("posterior cutoff:", toString(gam)))
  #gam <- gam + 0.005
  mem_results <- lapply(1:n_prior, function(prior_choice) {
    mem_fwer <- NULL
    for(j in 1:nrow(pmat)){
      pp <- mem_pp[prior_choice, j, , ]
      mem_pow[prior_choice, j, ] <- colMeans(pp > (gam[prior_choice]))
      mem_fwer[j] <- fwer1(post_prob = matrix(pp[, pmat[j, ]==p0], nrow=nsim),
                           gam = gam[prior_choice])$FWER
    }
    results <- cbind(mem_pow[prior_choice, , ], mem_fwer)
    colnames(results)[1:B] <- LETTERS[1:B]
    return(results)
    # print(cbind(mem_pow[prior_choice, , ], mem_fwer))
  })
  return(mem_results)
}