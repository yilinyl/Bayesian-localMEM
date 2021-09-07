# ---- Tuning stopping bounds ----

# pp_list <- lapply(1:10, function(i) {sim3_all[[1]]$pp_all[[i]]})

# qf <- 1 - 0.009^(n_enroll / Nmax_v)

# int_n <- c(10, 19)
# Nmax <- rep(19, 6)
# judgment for a single simulation
single_judge <- function(pp_mat, qf, n_stg, n_min = 10, n_mat = NULL) {
  n_basket <- nrow(pp_mat)
  is_futile <- setNames(rep(0, n_basket), LETTERS[1:n_basket])
  early_stop <- setNames(rep(0, n_basket), LETTERS[1:n_basket])
  
  # qf <- 1 - delta^int_p
  for (i in 1:ncol(pp_mat)) {
    id_tmp <- which(is_futile == 0)
    if (!is.null(n_mat)) {
      qf_tmp <- qf[n_mat[i, ]+1]
      # print(n_mat[i, ])
    }
    else {
      qf_tmp = rep(qf[i], n_basket)
    }
    is_futile[id_tmp][pp_mat[id_tmp, i] < qf_tmp[id_tmp]] <- 1
    if (i < n_stg) {  # early stopping for futility
      early_stop[id_tmp][pp_mat[id_tmp, i] < qf_tmp[id_tmp]] <- 1
    }
  }
  return(list("is_futile"=is_futile, "early_stop"=early_stop))
}

# pp_list: posterior probability of p > p0 (in every single simulation)
tuning <- function(pp_list, int_p=NULL, pvec, p0=0.15, n_list = NULL,  # modified
                   # lambda=0.95, 
                   n_min = 10, n_max = NULL, r_grid=seq(1, 2, 0.1), 
                   par_grid = NULL, qf_all = NULL,
                   d_grid = seq(0, 0.1, 0.001), fwer_thres = 0.1, 
                   single_stg = FALSE) {
  nsims <- length(pp_list)
  truth <- as.numeric(pvec > p0)
  if (is.null(int_p) & is.null(n_list)) {
    stop("interim size (n_list) should be specified when int_p is missing")
  }
  
  if (!is.null(n_list) & is.null(n_max)) {
    n_max <- max(sapply(n_list, max))
  }
  if (is.null(int_p)) {
    n_stg <- n_max - n_min + 1
    int_p <- seq(0, n_max, 1) / n_max
  }
  else {
    n_stg <- length(int_p)
  }
  
  specify_qf = FALSE
  pars <- NULL
  # if (length(int_p) > 1) {
  #   par_grid <- r_grid
  #   qf_all <- t(sapply(r_grid, function(r) {lambda * int_p^r}))
  # }
  # else {
  #   par_grid <- d_grid
  #   qf_all <- matrix(sapply(d_grid, function(d) {1 - d^int_p}))
  # }
  if (is.null(qf_all)) {
    if (is.null(par_grid)) {
      if (single_stg) {
        r_grid <- 0
        # par_grid <- matrix(d_grid, ncol = 1)
      }
      par_grid <- expand.grid(r_grid, d_grid)  # TODO: test BOP2 bound
    }
    
    print(paste("******Test", nrow(par_grid), "combinations******"), quote = F)
    qf_all <- t(sapply(1:nrow(par_grid), function(i)
      {(1-par_grid[i, 2]) * int_p^par_grid[i, 1]}))
  }
  else {
    specify_qf = TRUE
    print(paste("******Search for Qf directly [total", nrow(qf_all), "]******"), quote = F)
  }
  
  # par_grid <- d_grid
  # qf_all <- matrix(sapply(d_grid, function(d) {1 - d^int_p}),
  #                  ncol = length(int_p), byrow = T)
  pet_all <- NULL
  qf_record <- NULL
  for (i in 1:nrow(qf_all)) {
    qf <- qf_all[i, ]
    # tuning_par <- par_grid[i, ]  # TODO
    if (specify_qf) {
      tuning_par <- rep(NA, 2)
    }
    else {
      tuning_par <- unlist(par_grid[i, ])
    }
    
    judge_all <- lapply(1:length(pp_list), function(i) {
      n_mat = NULL
      if (!is.null(n_list)) {
        n_mat = n_list[[i]]
      }
      single_judge(pp_list[[i]], qf, n_stg, n_mat = n_mat)
      }) 
    outcomes <- sapply(judge_all, function(l) {l$is_futile}) # is_futile for all
    pet <- rowSums(sapply(judge_all, function(l) {l$early_stop})) / nsims  #PET
    
    rej_h0 <- rowSums(1 - outcomes) / nsims
    fwer <- NA
    # fwer_or_power <- NA
    if (all(truth == 0)) {  # global null: compute FWER
      fwer <- sum(apply(outcomes, 2, function(c) {as.numeric(any(c == 0))})) / nsims
      print(paste("Qf:", paste(round(qf, 3), collapse = ", "), "FWER:", round(fwer, 3)))
      # print(paste("delta:", tuning_par, "fwer:", fwer))
      if (fwer < fwer_thres) {
        qf_record <- rbind(qf_record, qf)
        pars <- rbind(pars, c(tuning_par, rej_h0, fwer))
        pet_all <- rbind(pet_all, c(tuning_par, pet))
      }
    }
    else {  # evaluate using fixed delta
      # power
      truth_map = as.logical(truth)
      # power1 <- sum(apply(outcomes, 2, function(record) {
      #   all(record[truth_map] == 0)
      # })) / nsims
      power2 <- sum(apply(outcomes, 2, function(record) {
        all(record[truth_map] == 0)
      })) / nsims
      # print()
      qf_record <- rbind(qf_record, qf)
      pars <- rbind(pars, c(tuning_par, rej_h0, power2))
      pet_all <- rbind(pet_all, c(tuning_par, pet))
    }
  }  # end for-loop
  colnames(pars) <- c(paste0("tuning_par", 1:length(tuning_par)),
                      paste0("power", LETTERS[1:6]), "fwer_or_power")
  return(list("pars"=pars, "pet"=pet_all, "qf"=qf_record))
  
}
# tune1 <- tuning(sim1_all[[1]]$pp_all, int_p = 1, pvec = pmat[1, ], fwer_thres = 0.1)
# tune3 <- tuning(sim3_all[[1]]$pp_all, int_p = c(10, 19)/19, pmat[1, ], fwer_thres = 0.1)
