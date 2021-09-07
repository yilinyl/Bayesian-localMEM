
# ---- interim analysis ----
# return: pp = posterior probability that p > p0
interim_analysis <- function(x_vec, n_vec, p0, bf_thres=3.2,...) {
  n_basket <- length(n_vec)
  if (n_basket == 1) {
    return(list("pp"=1 - pbeta(p0, 1 + x_vec, 1 + n_vec - x_vec),
                "is_pool"=NA))
  }
  partitions <- basket_part(x=x_vec, n=n_vec, max_cl = n_basket)
  # print(paste("BF threshold:", bf_thres))
  test <- basket_cluster(partitions, x=x_vec, n=n_vec, bf_thres = bf_thres)
  return(list("pp"=1 - pbeta(p0, test$alpha, test$beta),
              "is_pool"=test$is_pool))
}

# ---- monitoring a basket trial ----
#' @param props proportion of each subtype [vector]
#' @param Nmin minimum sample size for a basket to enter interim analysis
#' @param Nmax maximum sample size for each basket to finish interim analysis [vector]
#' @param step an interim look is conducted each time when additional `step` patients are enrolled
#' @param p0 response rate under H0
#' @param delta tuning parameter for cutoff values

# Feb.26 - same props for all, record pooling rate
sequential_monitor <- function(pvec, props, Nmin, Nmax_v, p0=0.15, step=5,
                               delta=0.02, r = 1, b_names=NULL, bmin=1, evaluate = TRUE,
                               verbose = FALSE, stepwise = FALSE, ...) {
  # Initialization
  n_basket <- length(props)
  if (is.null(b_names)) {
    b_names <- LETTERS[1:n_basket]
  }
  if (length(Nmin) == 1) {
    Nmin <- rep(Nmin, n_basket)
  }
  Nstep <- Nmin
  # initial recruitment
  enroll <- t(rmultinom(1, size = 1, prob = props))
  colnames(enroll) <- b_names
  n_enter <- sum(colSums(enroll) >= Nmin)
  prop_tmp <- props
  while (n_enter < bmin) {
    if (stepwise) {  # wait until all baskets reach Nmin
      prop_tmp[colSums(enroll) >= Nmin] <- 0
    }
    if (any(colSums(enroll) >= Nmax_v)) {
      props[colSums(enroll) >= Nmax_v] <- 0
      prop_tmp <- props
    }
    enroll <- rbind(enroll, t(rmultinom(1, size = 1, prob = prop_tmp)))
    n_enter <- sum(colSums(enroll) >= Nmin)
  }
  n_enroll <- colSums(enroll)
  # initial patient response data
  patient_dat <- lapply(1:length(pvec), function(i) {
    if (n_enroll[i] > 0) { rbinom(n_enroll[i], 1, pvec[i]) }
    else {0} })
  # --------------------------------------------
  is_futile <- setNames(rep(0, n_basket), b_names)
  continue <- setNames(rep(TRUE, n_basket), b_names)
  bounds <- data.frame()  # [qf | n | group]
  pool_record <- list("n_basket"=NULL, "is_pool"=NULL)
  x_all <- 0
  pp_record <- NULL
  n_anal <- 0
  x_record <- NULL
  n_record <- NULL
  # ---- outer loop: stop if all trials terminated (futile or exceed Nmax) ----
  while (sum(continue) > 0) {
    n_enroll <- colSums(enroll)
    # enter_id <- as.vector(which(n_enroll[continue] >= Nstep[continue]))
    # enter_id <- as.vector(which(continue)[(n_enroll>=Nstep)[continue]])
    is_enter <- continue & (n_enroll >= Nstep)
    b_remain <- b_names[continue]
    # n_vec <- n_enroll[continue][enter_id]
    n_vec <- n_enroll[is_enter]
    x_all <- x_all + sapply(1:n_basket, function(b) {sum(patient_dat[[b]])})
    x_vec <- x_all[is_enter]
    x_record <- rbind(x_record, x_all) 
    n_record <- rbind(n_record, n_enroll)
    # print("------------------------------------", quote = F)
    # print(paste("n_sample:", paste(n_vec, collapse = ",")))
    # print(paste("n_response:", paste(x_vec, collapse = ",")))
    # print("------------------------------------", quote = F)
    
    anal <- interim_analysis(x_vec, n_vec, p0)
    
    n_anal <- n_anal + 1
    pp_f <- anal$pp
    pool_record$n_basket <- c(pool_record$n_basket, length(n_vec))
    pool_record$is_pool <- c(pool_record$is_pool, anal$is_pool)
    if (n_anal == 1) {  # first analysis
      pp_tmp <- rep(1, n_basket)
    }
    else {
      pp_tmp <- pp_record[, n_anal-1]  # copy values from last interim look
    }
    # pp_tmp[b_names %in% b_remain] <- pp_f
    pp_tmp[is_enter] <- pp_f  # modified for accrual rate
    pp_record <- cbind(pp_record, pp_tmp)
    colnames(pp_record)[n_anal] <- n_anal
    
    if (!evaluate) {  # compute posterior probability only (no evaluation)
      delta <- 1
      # qf <- ((1 - delta)*(n_enroll / Nmax_v)^r)[continue][enter_id]
      # bounds <- rbind(bounds, cbind(as.vector(qf), as.vector(n_vec), names(qf)))
    }
    qf <- ((1 - delta)*(n_enroll / Nmax_v)^r)[is_enter]
    bounds <- rbind(bounds, cbind(as.vector(qf), as.vector(n_vec), names(qf)))
    
    # pp_mat[j, b_remain] <- as.vector(pp_f)
    # pp_mat[j, !continue] <- 1
    # pp_f, qf sliced only for continue arms that entered current analysis
    if (any(pp_f < qf)) {
      if (verbose == TRUE) {
        print(paste("* stop basket", paste(b_names[is_enter][pp_f < qf], collapse = ","), 
                    "for futility"))
      }
      is_futile[b_names[is_enter][pp_f < qf]] <- 1  # pp_f < delta: stop for f
      continue[b_names[is_enter][pp_f < qf]] <- FALSE
    }
    if (any(n_enroll >= Nmax_v)) {
      b_complete <- b_names[n_enroll >= Nmax_v]
      continue[b_names[n_enroll >= Nmax_v]] <- FALSE
      if (verbose & any(b_remain %in% b_complete)) {
        print(paste("* basket", paste(b_remain[b_remain %in% b_complete], collapse = ","), 
                    "reached Nmax <complete!>"))
      }
      
    }
    # if (sum(continue) == 0) {  # all trials are terminated
    #   break
    # }
    # continue recruitment
    props[!continue] <- 0  # stop recruitment for terminated baskets
    # Nstep <- Nmin
    if (stepwise) {
      Nstep <- Nstep + step
    }
    while (any(continue)) {  
      for (i in 1:step) {
        if (sum(props) == 0) {
          break
        }
        enroll <- rbind(enroll, t(rmultinom(1, 1, prob = props)))
        if (any(colSums(enroll) >= Nmax_v)) {
          props[colSums(enroll) >= Nmax_v] <- 0  # stop recruitment for baskets reaching max sample size
        }
      }
      # continue recruitment only if no basket could enter interim analysis
      if (sum(continue) > 1) {
        if (stepwise) {
          flag <- all(colSums(enroll[, continue]) >= Nstep[continue])
        }
        else {
          enroll_tmp <- colSums(enroll)
          enter_tmp <- as.vector(which(continue)[(enroll_tmp >= Nmin)[continue]])
          # flag <- any(colSums(enroll[, continue]) >= Nmin[continue])
          flag <- (length(enter_tmp) > 0) & 
            any(enroll_tmp[enter_tmp] != n_enroll[enter_tmp])  # avoid duplicate assessment
        }
      }
      else {
        flag <- sum(enroll[, continue]) >= Nmin[continue]
      }
      # flag <- ifelse(sum(continue) > 1,  
      #                # any(colSums(enroll[, continue]) > Nmin), 
      #                all(colSums(enroll[, continue]) > Nmin), 
      #                sum(enroll[, continue]) > Nmin)
      if (flag) { break }
    }
    n_new <- colSums(enroll) - n_enroll
    patient_dat <- lapply(1:length(pvec), function(i) {
      if (n_new[i] > 0) { rbinom(n_new[i], 1, pvec[i]) }
      else {0} 
      })
  }
  names(bounds) <- c("Qf", "n", "group")
  bounds$Qf <- as.numeric(bounds$Qf)
  bounds$n <- as.numeric(bounds$n)
  bounds <- unique(bounds)
  bounds_sep <- split(bounds, bounds$group)
  return(list("ptrue" = pvec, "is_futile" = is_futile, "bounds" = bounds_sep, 
              "responses"=x_record, "interim_n"=n_record,
              "n_enroll" = n_enroll, "pool_record" = data.frame(pool_record),
              "pp" = pp_record))
}

sim_and_eval <- function(pvec, Nvec, p0=0.15, delta = 0.2, Nmin = NULL,
                   props = NULL, nsims=50, bmin=1, step=5, stepwise = FALSE,
                   seed=123, print_freq=10,...) {
  n_basket <- length(pvec)
  truth <- as.numeric(pvec > p0)
  set.seed(seed)
  if (is.null(props)) {
    props <- Nvec / sum(Nvec)
  }
  if (is.null(Nmin)) {
    Nmin <- min(Nvec)
  }
  out <- lapply(1:nsims, function(s) {
    if (s %% print_freq == 0) {
      print(paste("---- Simulation", s, "/", nsims, "----"))
    }
    sequential_monitor(pvec, props, Nmin, Nmax_v = Nvec, 
                       delta = delta, bmin = bmin, step = step, 
                       stepwise = stepwise)
  })
  outcomes <- sapply(out, function(l) {l$is_futile})
  pp_record <- lapply(out, function(l) {l$pp})  # record pp
  rej_h0 <- rowSums(1 - outcomes) / nsims
  fwer <- NA
  if (all(truth == 0)) {  # global null: compute FWER
    fwer <- sum(apply(outcomes, 2, function(c) {as.numeric(any(c == 0))})) / nsims
  }
  n_looks <- sum(sapply(out, function(l) {nrow(na.omit(l$pool_record))}))
  p_pool <- sum(sapply(out, function(l) 
    {sum(na.omit(l$pool_record)$is_pool)})) / n_looks
  
  return(list("pvec"=pvec, "marginal_power" = rej_h0, "p0"=p0, "FWER"=fwer,
              "p_pool" = p_pool, "delta" = delta, "pp_all"=pp_record))
  # "Qf" = 1 - delta^(Nmin:Nvec[1]/Nvec[1])
}

# int_n <- c(10, 9)
pp_direct <- function(pvec, int_n, p0=0.15, b_names=NULL) {
  n_basket <- length(pvec)
  # Nmax <- sum(int_n)
  if (is.null(b_names)) {
    b_names <- LETTERS[1:n_basket]
  }
  if (is.null(dim(int_n))) {
    int_n <- matrix(rep(int_n, n_basket), ncol = n_basket)
  }
  xmat <- t(apply(int_n, 1, function(nvec) {rbinom(n_basket, nvec, pvec)}))
  Nmat <- apply(int_n, 2, cumsum)
  out <- lapply(1:nrow(xmat), function(i) 
    {interim_analysis(xmat[i, ], Nmat[i, ], p0=0.15)})
  pp <- sapply(out, function(l) {l$pp})
  pool_record <- sapply(out, function(l) {l$is_pool})
  return(list("pp" = pp, "is_pool" = pool_record))
}

# lapply(1:10, function(s) {pp_direct(pvec, c(10, 9))})
