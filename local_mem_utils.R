# local MEM utils
require(partitions)

calc_priors <- function(n_basket, n_priors, max_cl=n_basket, part=NULL) {
  if (is.null(part)) {
    part <- get.part(R = n_basket, max_cl = n_basket)
  }
  K <- nrow(part)
  n_bk <- apply(part, 1, function(x){
    length(unique(x))
  })
  
  prior_mat <- matrix(NA, nrow = K, ncol = n_priors)
  # set the prior choice
  for (prior_choice in 1:n_priors) {
    prior_mat[, prior_choice] <- n_bk**(prior_choice-1) / sum(n_bk**(prior_choice-1))
  }
  # print(colSums(prior_mat))
  return(cbind(part, prior_mat))
}

# calc_priors(4, 3)
