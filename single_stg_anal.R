# single stage analysis

single_stg_anal <- function(pvec, pp, qf_grid, p0=0.15, alpha=0.1, tuning=FALSE) {
  # truth <- as.numeric(pvec > p0)
  n_basket <- ncol(pp)
  nsims <- nrow(pp)
  pars <- NULL
  for (qf in qf_grid) {
    is_futile <- apply(pp, 1, function(pp_i) {as.numeric(pp_i < qf)})
    rej_h0 <- rowSums(1 - is_futile) / nsims  # typeI error(H0) / marginal power(H1)
    # fwer <- NA
    
    # fwer <- sum(sapply(1:ncol(is_futile), function(i) {
    #   any(is_futile[, i] == 0)
    # })) / nsims
    fwer <- mean(apply(is_futile, 2, function(is_fut) {
      any(is_fut[pvec <= p0]==0)}))
    if (tuning) {
      if (fwer < alpha) {
        pars <- rbind(pars, c(qf, rej_h0, fwer))
      }
    }
    else {
      pars <- rbind(pars, c(qf, rej_h0, fwer))
    }
    
  }  # end for-loop
  colnames(pars) <- c("qf", LETTERS[1:n_basket], "fwer")
  return(pars)
}