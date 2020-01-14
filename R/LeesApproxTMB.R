
logit <- function(x) qlogis(x)

interpolation_in_interval <- function(v, x, cond) {
  if(cond == 1) {
    # If v == x Returns the index for yout (remember C++ indexing starts at zero),
    # else returns a negative number if v is outside the range of x (LAA)
    if(sum(v == x) == 1) {
      return(which(v == x)-1)
    } else if(v < min(x) || v > max(x)) {
      return(-1)
    } else return(-2)
  }

  if(cond == 2 && all(v != x) && v > min(x) && v < max(x)) return(max(which(v > x)) - 1) # x1 in TMB code, x2 is x1 + 1

  if(cond == 3 && all(v != x) && v > min(x) && v < max(x)) {
    test <- c(max(which(v > x)), min(which(v < x)))
    return(test[2] - test[1] == 1)
  }

  return(NA)
}

find_interval <- function(xout, LenBins, Nbins) xout >= LenBins[1:Nbins] & xout <= LenBins[2:(Nbins+1)]

make_integrate_indices <- function(ind) {
  res <- apply(ind, c(1, 3), function(x) which(as.logical(x))-1)
  res2 <- lapply(res, identity)

  f_a_l <-  rep(1:length(res), times = vapply(res2, length, numeric(1))) - 1
  integ_ind <- do.call(c, res)
  return(list(f_a_l, integ_ind))
}

#' Lees Approximation interpolation function
#'
#' Primarily used to verify that output from TMB will match that from Rcpp to numerical precision.
#'
#' @import TMB
#' @useDynLib LeesApproxTMB
#' @export
LeesApproxTMB <- function(FVec, ngtg, maxsd, binwidth, M,
                          Linf, K, t0, LFS, L5, Vmaxlen, LinfCV, maxage) {

  if(length(FVec) == 1) FVec <- rep(FVec, maxage)
  if(length(FVec) > maxage) FVec <- FVec[(length(FVec)-maxage+1):length(FVec)]

  distGTG <- seq(-maxsd, maxsd, length.out = ngtg)
  rdist <- dnorm(distGTG)
  rdist <- rdist/sum(rdist)
  Linfgtg <- Linf + LinfCV * Linf * distGTG

  LenBins <- seq(0, max(Linfgtg), binwidth)
  LenMids <- LenBins - 0.5 * binwidth
  LenMids <- LenMids[LenMids > 0]

  LAA <- matrix(0, maxage, ngtg)
  ages <- 1:maxage
  for(g in 1:ngtg) LAA[, g] <- Linfgtg[g] * (1 - exp(-K * (ages - t0)))

  xout <- t(apply(LAA, 1, function(x, y) sort(c(x, y)), y = LenBins))

  # For interpolation function int_point
  interp_check <- interp_check2 <- matrix(NA, maxage, ncol(xout))
  for(a in 1:maxage) {
    interp_check[a, ] <- vapply(xout[a, ], interpolation_in_interval, numeric(1), x = LAA[a, ], cond = 1)
    interp_check2[a, ] <- vapply(xout[a, ], interpolation_in_interval, numeric(1), x = LAA[a, ], cond = 2)
  }

  # For integration
  Nbins <- length(LenMids)
  ind <- array(NA, c(Nbins, ncol(xout), maxage))
  for(a in 1:maxage) {
    ind[, , a] <- vapply(xout[a, ], find_interval, numeric(Nbins), LenBins = LenBins, Nbins = Nbins)
  }
  integ_check <- apply(ind, c(3, 1), function(x) as.integer(sum(x) == 1))
  integ_ind <- make_integrate_indices(ind)

  TMB_data <- list(model = "LeesApprox_internal",
                   FVec = FVec, ngtg = ngtg, LenBins = LenBins, LenMids = LenMids, ages = ages,
                   M = M, Linf = Linf, LAA = LAA, xout = xout, LFS = LFS, L5 = L5, Vmaxlen = Vmaxlen,
                   maxage = maxage, distGTG = distGTG, rdist = rdist,
                   interp_check = interp_check, interp_check2 = interp_check2, integ_check = integ_check,
                   integ_fac = integ_ind[[1]], integ_ind = integ_ind[[2]])

  obj <- MakeADFun(data = TMB_data, parameters = list(p = 0), silent = TRUE, DLL = "LeesApproxTMB")
  report <- obj$report()

  return(list(probLA = report$probLA, LenComp = report$LenComp, Select_at_age = report$Select_at_age,
              Select_at_length = report$Select_at_length, LAA = report$LAA, LenMids = report$LenMids,
              LenBins = report$LenBins, NAA = report$NAA))
}

#' SCA with Growth Type Groups (GTG)
#'
#' Includes interpolation approximation for Lee's effect.
#' @param x A position in the Data object (by default, equal to one for assessments).
#' @param Data An object of class Data
#' @param ngtg Integer, the number of growth type groups in the model.
#' @param max_sd_gtg The number of standard deviations to span the GTGs across the length-at-age distribution.
#' @param use_LeesEffect Logical, whether to incorporate Lee's effect in the model.
#' @param truncate_CAL Logical, whether to truncate the length composition data at the tails if the length bins are outside the
#' range of lengths of all growth type groups.
#' @param CAL_multiplier Numeric for data weighting of catch-at-length matrix.
#' Default is set to zero to ignore length comp data.
#' @param CAA_multiplier Numeric for data weighting of catch-at-age matrix.
#' Default is set to zero to ignore age comp data.
#' @param SR Stock-recruit function (either \code{"BH"} for Beverton-Holt or \code{"Ricker"}).
#' @param vulnerability Whether estimated vulnerability is \code{"logistic"} or \code{"dome"} (double-normal).
#' See details for parameterization.
#' @param I_type Whether the index surveys population biomass (B; this is the default in the DLMtool operating model),
#' vulnerable biomass (VB), or spawning stock biomass (SSB).
#' @param rescale A multiplicative factor that rescales the catch in the assessment model, which
#' can improve convergence. By default, \code{"mean1"} scales the catch so that time series mean is 1, otherwise a numeric.
#' Output is re-converted back to original units.
#' @param max_age Integer, the maximum age (plus-group) in the model.
#' @param start Optional list of starting values. Entries can be expressions that are evaluated in the function. See details.
#' @param fix_h Logical, whether to fix steepness to value in \code{Data@@steep} in the model for \code{SCA}. This only affects
#' calculation of reference points for \code{SCA2}.
#' @param fix_F_equilibrium Logical, whether the equilibrium fishing mortality prior to the first year of the model
#' is estimated. If \code{TRUE}, \code{F_equilibrium} is fixed to value provided in \code{start} (if provided),
#' otherwise, equal to zero (assumes unfished conditions).
#' @param fix_U_equilibrium Logical, same as `fix_F_equilibrium` for `SCA_Pope`.
#' @param fix_omega Logical, whether the standard deviation of the catch is fixed. If \code{TRUE},
#' sigma is fixed to value provided in \code{start} (if provided), otherwise, value based on \code{Data@@CV_Cat}.
#' @param fix_sigma Logical, whether the standard deviation of the index is fixed. If \code{TRUE},
#' sigma is fixed to value provided in \code{start} (if provided), otherwise, value based on \code{Data@@CV_Ind}.
#' @param fix_tau Logical, the standard deviation of the recruitment deviations is fixed. If \code{TRUE},
#' tau is fixed to value provided in \code{start} (if provided), otherwise, value based on \code{Data@@sigmaR}.
#' @param common_dev Typically, a numeric for the number of most recent years in which a common recruitment deviation will
#' be estimated (in \code{SCA2}, uninformative years will have a recruitment closer to the mean, which can be very misleading,
#' especially near the end of the time series). By default, \code{"comp50"} uses the number of ages (smaller than the mode)
#' for which the catch-at-age matrix has less than half the abundance than that at the mode.
#' @param early_dev Character string describing the years for which recruitment deviations are estimated in \code{SCA}. By default, \code{"comp_onegen"}
#' rec devs are estimated one full generation prior to the first year when catch-at-age (CAA) data are available. With \code{"comp"}, rec devs are
#' estimated starting in the first year with CAA. With \code{"all"}, rec devs start at the beginning of the model.
#' @param late_dev Typically, a numeric for the number of most recent years in which recruitment deviations will
#' not be estimated in \code{SCA} (recruitment in these years will be based on the mean predicted by stock-recruit relationship).
#' By default, \code{"comp50"} uses the number of ages (smaller than the mode)
#' for which the catch-at-age matrix has less than half the abundance than that at the mode.
#' @param integrate Logical, whether the likelihood of the model integrates over the likelihood
#' of the recruitment deviations (thus, treating it as a random effects/state-space variable).
#' Otherwise, recruitment deviations are penalized parameters.
#' @param silent Logical, passed to \code{\link[TMB]{MakeADFun}}, whether TMB
#' will print trace information during optimization. Used for dignostics for model convergence.
#' @param opt_hess Logical, whether the hessian function will be passed to \code{\link[stats]{nlminb}} during optimization
#' (this generally reduces the number of iterations to convergence, but is memory and time intensive and does not guarantee an increase
#' in convergence rate). Ignored if \code{integrate = TRUE}.
#' @param n_restart The number of restarts (calls to \code{\link[stats]{nlminb}}) in the optimization procedure, so long as the model
#' hasn't converged. The optimization continues from the parameters from the previous (re)start.
#' @param control A named list of agruments for optimization to be passed to
#' \code{\link[stats]{nlminb}}.
#' @param inner.control A named list of arguments for optimization of the random effects, which
#' is passed on to \code{\link[TMB]{newton}}.
#' @param ... Other arguments to be passed.
#' @import TMB MSEtool
#' @useDynLib LeesApproxTMB
#' @export
SCA_GTG <- function(x = 1, Data, SR = c("BH", "Ricker"), vulnerability = c("logistic", "dome"),
                    CAA_multiplier = 0, CAL_multiplier = 25, ngtg = 3L, max_sd_gtg = 2, use_LeesEffect = TRUE, truncate_CAL = TRUE,
                    I_type = c("B", "VB", "SSB"), rescale = "mean1", max_age = Data@MaxAge,
                    start = NULL, fix_h = TRUE, fix_F_equilibrium = TRUE, fix_omega = TRUE, fix_sigma = FALSE, fix_tau = TRUE,
                    early_dev = c("comp_onegen", "comp", "all"), late_dev = "comp50", integrate = FALSE,
                    silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                    control = list(iter.max = 2e5, eval.max = 4e5), inner.control = list(), ...) {

  dependencies <- "Data@Cat, Data@Ind, Data@Mort, Data@L50, Data@L95, Data@CAL, Data@vbK, Data@vbLinf, Data@LenCV, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge"
  dots <- list(...)
  start <- lapply(start, eval, envir = environment())

  max_age <- as.integer(min(max_age, Data@MaxAge))
  vulnerability <- match.arg(vulnerability)
  SR <- match.arg(SR)
  I_type <- match.arg(I_type)
  early_dev <- match.arg(early_dev)

  if(any(names(dots) == "yind")) {
    yind <- eval(dots$yind)
  } else {
    yind <- which(!is.na(Data@Cat[x, ]))[1]
    yind <- yind:length(Data@Cat[x, ])
  }
  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]
  if(any(is.na(C_hist) | C_hist < 0)) warning("Error. Catch time series is not complete.")
  I_hist <- Data@Ind[x, yind]

  n_y <- length(C_hist)
  M <- rep(Data@Mort[x], max_age)
  a <- Data@wla[x]
  b <- Data@wlb[x]
  Linf <- Data@vbLinf[x]
  K <- Data@vbK[x]
  t0 <- Data@vbt0[x]
  La <- Linf * (1 - exp(-K * (c(1:max_age) - t0)))
  Wa <- a * La ^ b
  A50 <- min(0.5 * max_age, iVB(t0, K, Linf, Data@L50[x]))
  A95 <- max(A50+0.5, iVB(t0, K, Linf, Data@L95[x]))
  mat_age <- 1/(1 + exp(-log(19) * (c(1:max_age) - A50)/(A95 - A50)))
  mat_age <- mat_age/max(mat_age)

  # GTG
  distGTG <- seq(-max_sd_gtg, max_sd_gtg, length.out = ngtg)
  rdist <- dnorm(distGTG)
  rdist <- rdist/sum(rdist)
  LenCV <- Data@LenCV[x]
  Linfgtg <- Linf + LenCV * Linf * distGTG

  LAA <- outer(1 - exp(-K * (c(1:max_age) - t0)), Linfgtg)

  # Age comps
  CAA_hist <- try(Data@CAA[x, yind, 1:max_age], silent = TRUE)
  if(is.character(CAA_hist)) CAA_hist <- matrix(NA, length(yind), max_age)
  if(max_age < Data@MaxAge) CAA_hist[, max_age] <- rowSums(Data@CAA[x, yind, max_age:Data@MaxAge], na.rm = TRUE)

  CAA_n_nominal <- rowSums(CAA_hist)
  if(CAA_multiplier <= 1) {
    CAA_n_rescale <- CAA_multiplier * CAA_n_nominal
  } else CAA_n_rescale <- pmin(CAA_multiplier, CAA_n_nominal)
  CAA_n_rescale[!CAA_n_rescale] <- NA

  # Length comps
  CAL_hist <- try(Data@CAL[x, yind, ], silent = TRUE)
  if(is.character(CAL_hist)) CAL_hist <- matrix(NA, length(yind), length(Data@CAL_mids))
  CAL_bins <- Data@CAL_bins
  CAL_mids <- Data@CAL_mids

  if(truncate_CAL) {
    CAL_min <- ifelse(all(!Data@CAL_bins < min(LAA)), 1, max(which(Data@CAL_bins < min(LAA))))
    CAL_max <- ifelse(all(!Data@CAL_bins > max(LAA)), length(Data@CAL_bins), min(which(Data@CAL_bins > max(LAA))))

    CAL_hist <- CAL_hist[ , CAL_min:(CAL_max-1), drop = FALSE]
    CAL_bins <- CAL_bins[CAL_min:CAL_max]
    CAL_mids <- CAL_mids[CAL_min:(CAL_max-1)]
  }

  CAL_n_nominal <- rowSums(CAL_hist)
  if(CAL_multiplier <= 1) {
    CAL_n_rescale <- CAL_multiplier * CAL_n_nominal
  } else CAL_n_rescale <- pmin(CAL_multiplier, CAL_n_nominal)

  xout <- t(apply(LAA, 1, function(x, y) sort(c(x, y)), y = CAL_bins))

  # Interpolation
  # Determine the GTG indices for interpolated abundance at each xout (ordered concatenation of LenBin and LAA for each GTG)
  # interp_check evaluates if the j-th entry of xout = one of LAA or outside the range of LAA, otherwise
  # interp_check2 identifies the indices of LAA in which xout is in between
  interp_check <- interp_check2 <- matrix(NA, max_age, ncol(xout))
  for(aa in 1:max_age) {
    interp_check[aa, ] <- vapply(xout[aa, ], interpolation_in_interval, numeric(1), x = LAA[aa, ], cond = 1)
    interp_check2[aa, ] <- vapply(xout[aa, ], interpolation_in_interval, numeric(1), x = LAA[aa, ], cond = 2)
  }

  # For integration of abundance for each length bin
  # ind identifies the indices of xout within each length bin
  # integ_check evaluates whether xout is on the boundary of the smallest/largest length bin, otherwise
  # integ_ind returns the indices of yout
  # integ_ind is a list with the indices corresponding to age and length and the indices for each xout/yout
  # over which to sum abundances
  Nbins <- length(CAL_mids)
  ind <- array(NA, c(Nbins, ncol(xout), max_age))
  for(aa in 1:max_age) {
    ind[, , aa] <- vapply(xout[aa, ], find_interval, numeric(Nbins), LenBins = CAL_bins, Nbins = Nbins)
  }
  integ_check <- apply(ind, c(3, 1), function(x) as.integer(sum(x) == 1))
  integ_ind <- make_integrate_indices(ind)

  LH <- list(mean_LAA = La, mean_WAA = Wa, Linf = Linf, K = K, t0 = t0, a = a, b = b, A50 = A50, A95 = A95, LenCV = LenCV)

  if(early_dev == "all") {
    est_early_rec_dev <- rep(1, max_age-1)
    est_rec_dev <- rep(1, n_y)
  }
  if(early_dev == "comp") {
    est_early_rec_dev <- rep(NA, max_age-1)
    ind1 <- which(!is.na(CAL_n_nominal))[1]
    est_rec_dev <- ifelse(c(1:n_y) < ind1, NA, 1)
  }
  if(early_dev == "comp_onegen") {
    ind1 <- which(!is.na(CAL_n_nominal))[1] - max_age
    if(ind1 < 0) {
      early_start <- max_age + ind1
      est_early_rec_dev <- rev(ifelse(c(1:(max_age-1)) < early_start, NA, 1))
      est_rec_dev <- rep(1, n_y)
    } else {
      est_early_rec_dev <- rep(NA, max_age-1)
      est_rec_dev <- ifelse(c(1:n_y) < ind1, NA, 1)
    }
  }
  if(is.character(late_dev) && late_dev == "comp50") {
    CAL_all <- colSums(CAL_hist, na.rm = TRUE)/max(colSums(CAL_hist, na.rm = TRUE))
    CAL_mode <- which.max(CAL_all)[1]
    comp50_ind <- which(CAL_all[1:CAL_mode] <= 0.5)
    comp50_ind <- comp50_ind[length(comp50_ind)]
    late_dev <- ifelse(is.na(comp50_ind), 0, which.min(abs(CAL_mids[comp50_ind] - La)))
  }
  if(is.numeric(late_dev) && late_dev > 0) {
    if(late_dev > length(est_rec_dev)) late_dev <- length(est_rec_dev)
    ind_late <- (length(est_rec_dev) - late_dev + 1):length(est_rec_dev)
    est_rec_dev[ind_late] <- NA
  }

  if(rescale == "mean1") rescale <- 1/mean(C_hist)
  data <- list(model = "LeesApprox_SCA", C_hist = C_hist * rescale, I_hist = I_hist,
               CAA_hist = t(apply(CAA_hist, 1, function(x) x/sum(x))), CAA_n = CAA_n_rescale,
               CAL_hist = t(apply(CAL_hist, 1, function(x) x/sum(x))), CAL_n = CAL_n_rescale,
               n_y = n_y, max_age = max_age, M = M,
               WAA = a * LAA^b, mat = mat_age, I_type = I_type,
               SR_type = SR, est_early_rec_dev = est_early_rec_dev, est_rec_dev = est_rec_dev,
               ngtg = ngtg, Nbins = Nbins, LenBins = CAL_bins, LenMids = CAL_mids, Linf = Linf,
               LAA = LAA, xout = xout, distGTG = distGTG, rdist = rdist,
               interp_check = interp_check, interp_check2 = interp_check2, integ_check = integ_check,
               integ_fac = integ_ind[[1]], integ_ind = integ_ind[[2]],
               use_LeesEffect = as.integer(use_LeesEffect), yind_F = as.integer(0.5 * n_y))
  data$CAA_hist[data$CAA_hist < 1e-8] <- 1e-8
  data$CAL_hist[data$CAL_hist < 1e-8] <- 1e-8

  # Starting values
  params <- list()
  if(!is.null(start)) {
    if(!is.null(start$R0) && is.numeric(start$R0)) params$log_R0 <- log(start$R0[1] * rescale)
    if(!is.null(start$h) && is.numeric(start$h)) {
      if(SR == "BH") {
        h_start <- (start$h[1] - 0.2)/0.8
        params$transformed_h <- logit(h_start)
      } else {
        params$transformed_h <- log(start$h[1] - 0.2)
      }
    }
    if(!is.null(start$F_equilibrium) && is.numeric(start$F_equilibrium)) params$F_equilibrium <- start$F_equilibrium
    if(!is.null(start$vul_par) && is.numeric(start$vul_par)) {
      if(start$vul_par[1] > 0.9 * Linf) stop("start$vul_par[1] needs to be less than 0.9 * Data@Linf.")
      if(vulnerability == "logistic") {
        if(length(start$vul_par) < 2) stop("Two parameters needed for start$vul_par with logistic vulnerability (see help).")
        if(start$vul_par[1] <= start$vul_par[2]) stop("start$vul_par[1] needs to be greater than start$vul_par[2] (see help).")

        params$vul_par <- c(logit(start$vul_par[1]/Linf/0.9), log(start$vul_par[1] - start$vul_par[2]), 10)
      }
      if(vulnerability == "dome") {
        if(length(start$vul_par) < 3) stop("Three parameters needed for start$vul_par with dome vulnerability (see help).")
        if(start$vul_par[1] <= start$vul_par[2]) stop("start$vul_par[1] needs to be greater than start$vul_par[2] (see help).")
        if(start$vul_par[3] <= 0 || start$vul_par[3] >= 1) stop("start$vul_par[3] needs to be between 0-1 (see help).")

        params$vul_par <- c(logit(start$vul_par[1]/Linf/0.9), log(start$vul_par[1] - start$vul_par[2]),
                            logit(start$vul_par[3]))
      }
    }
    if(!is.null(start$F) && is.numeric(start$F)) params$logF <- log(start$F)

    if(!is.null(start$omega) && is.numeric(start$omega)) params$log_omega <- log(start$omega)
    if(!is.null(start$sigma) && is.numeric(start$sigma)) params$log_sigma <- log(start$sigma)
    if(!is.null(start$tau) && is.numeric(start$tau)) params$log_tau <- log(start$tau)
  }

  if(is.null(params$log_R0)) {
    params$log_R0 <- ifelse(is.null(Data@OM$N0[x]), log(mean(data$C_hist)) + 4,
                            log(1.5 * rescale * Data@OM$R0[x]))
  }
  if(is.null(params$transformed_h)) {
    h_start <- ifelse(!fix_h && is.na(Data@steep[x]), 0.9, Data@steep[x])
    if(SR == "BH") {
      h_start <- (h_start - 0.2)/0.8
      params$transformed_h <- logit(h_start)
    } else {
      params$transformed_h <- log(h_start - 0.2)
    }
  }
  if(is.null(params$F_equilibrium)) params$F_equilibrium <- 0
  if(is.null(params$vul_par)) {
    if((is.na(Data@LFC[x]) && is.na(Data@LFS[x])) || (Data@LFC[x] > Linf) || (Data@LFS[x] > Linf)) {
      CAL_mode <- which.max(colSums(CAL_hist, na.rm = TRUE))
      if(all(is.na(CAL_hist))) {
        CAL_mode <- La[which.max(colSums(CAA_hist, na.rm = TRUE))] 
      }
      LFS <- CAL_mids[CAL_mode]
      L5 <- 0.4 * LFS
    } else {
      LFS <- Data@LFS[x]
      L5 <- Data@LFC[x]
    }
    if(vulnerability == "logistic") params$vul_par <- c(logit(LFS/Linf/0.9), log(LFS - L5), 10)
    if(vulnerability == "dome") {
      params$vul_par <- c(logit(LFS/Linf/0.9), log(LFS - L5), logit(0.5))
    }
  }
  if(is.null(params$logF)) {
    logFstart <- numeric(n_y)
    logFstart[data$yind_F + 1] <- 0.75 * mean(data$M)
    params$logF <- logFstart
  }

  if(is.null(params$log_sigma)) {
    sigmaI <- max(0.01, sdconv(1, Data@CV_Ind[x]), na.rm = TRUE)
    params$log_sigma <- log(sigmaI)
  }
  if(is.null(params$log_omega)) {
    sigmaC <- max(0.01, sdconv(1, Data@CV_Cat[x]), na.rm = TRUE)
    params$log_omega <- log(sigmaC)
  }
  if(is.null(params$log_tau)) {
    tau_start <- ifelse(is.na(Data@sigmaR[x]), 0.6, Data@sigmaR[x])
    params$log_tau <- log(tau_start)
  }
  params$log_early_rec_dev <- rep(0, max_age - 1)
  params$log_rec_dev <- rep(0, n_y)

  info <- list(Year = Year, data = data, params = params, LH = LH, control = control,
               inner.control = inner.control, rescale = rescale)

  map <- list()
  if(any(info$data$C_hist <= 0)) {
    ind <- info$data$C_hist <= 0
    info$params$logF[ind] <- -20
    map_logF <- length(params$logF)
    map_logF[ind] <- NA
    map_logF[!ind] <- 1:sum(!ind)
    map$logF <- factor(map_logF)
  }
  if(fix_h) map$transformed_h <- factor(NA)
  if(fix_F_equilibrium) map$F_equilibrium <- factor(NA)
  if(fix_sigma) map$log_sigma <- factor(NA)
  if(fix_omega) map$log_omega <- factor(NA)
  if(fix_tau) map$log_tau <- factor(NA)
  if(vulnerability == "logistic") map$vul_par <- factor(c(1, 2, NA))
  if(any(is.na(est_early_rec_dev))) {
    n_est <- sum(!is.na(est_early_rec_dev))
    if(n_est == 0) map$log_early_rec_dev <- factor(rep(NA, max_age - 1))
    else {
      est_early_rec_dev[!is.na(est_early_rec_dev)] <- 1:n_est
      map$log_early_rec_dev <- factor(est_early_rec_dev)
    }
  }
  if(any(is.na(est_rec_dev))) {
    n_est <- sum(!is.na(est_rec_dev))
    est_rec_dev[!is.na(est_rec_dev)] <- 1:n_est
    map$log_rec_dev <- factor(est_rec_dev)
  }

  random <- NULL
  if(integrate) random <- c("log_early_rec_dev", "log_rec_dev")

  obj <- MakeADFun(data = info$data, parameters = info$params, hessian = TRUE,
                   map = map, random = random, DLL = "LeesApproxTMB", inner.control = inner.control, silent = silent)

  mod <- MSEtool:::optimize_TMB_model(obj, control, opt_hess, n_restart)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best)
  report$Len_at_age <- t(sapply(report$NPR, function(x, y) rowSums(x*y)/rowSums(x), y = info$data$LAA))
  report$Wt_at_age <- t(sapply(report$NPR, function(x, y) rowSums(x*y)/rowSums(x), y = info$data$WAA))
  report$VB0 <- NA

  if(rescale != 1) {
    vars_div <- c("B", "E", "CAApred", "CALpred", "CN", "Cpred", "N", "VB", "R", "R_early", "R_eq", "VB0", "R0", "B0", "E0", "N0")
    vars_mult <- "Brec"
    var_trans <- c("R0", "q")
    fun_trans <- c("/", "*")
    fun_fixed <- c("log", NA)
    MSEtool:::rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
  }

  Yearplusone <- c(Year, max(Year) + 1)
  YearEarly <- (Year[1] - max_age + 1):(Year[1] - 1)
  YearDev <- c(YearEarly, Year)
  YearR <- c(YearDev, max(YearDev) + 1)
  R <- c(rev(report$R_early), report$R)

  Dev <- c(rev(report$log_early_rec_dev), report$log_rec_dev)
  Dev_out <- structure(Dev, names = YearDev)

  nll_report <- ifelse(is.character(opt), ifelse(integrate, NA, report$nll), opt$objective)
  Assessment <- new("Assessment", Model = "SCA_GTG", Name = Data@Name, conv = !is.character(SD) && SD$pdHess,
                    B0 = report$B0, R0 = report$R0, N0 = report$N0,
                    SSB0 = report$E0, VB0 = report$VB0,
                    h = report$h, FMort = structure(report$F, names = Year),
                    B = structure(report$B, names = Yearplusone),
                    B_B0 = structure(report$B/report$B0, names = Yearplusone),
                    SSB = structure(report$E, names = Yearplusone),
                    SSB_SSB0 = structure(report$E/report$E0, names = Yearplusone),
                    VB = structure(report$VB, names = Yearplusone),
                    #VB_VB0 = structure(report$VB/report$VB0, names = Yearplusone),
                    R = structure(R, names = YearR),
                    N = structure(rowSums(report$N), names = Yearplusone),
                    N_at_age = report$N,
                    Selectivity = report$Select_at_age,
                    Obs_Catch = structure(C_hist, names = Year),
                    Obs_Index = structure(I_hist, names = Year),
                    Obs_C_at_age = CAA_hist,
                    Catch = structure(report$Cpred, names = Year),
                    Index = structure(report$Ipred, names = Year),
                    C_at_age = report$CAApred,
                    Dev = Dev_out,
                    Dev_type = "log-Recruitment deviations",
                    NLL = structure(c(nll_report, report$nll_comp, report$penalty + report$prior),
                                    names = c("Total", "Index", "CAA", "CAL", "Catch", "Dev", "Penalty")),
                    info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                    dependencies = dependencies)

  if(Assessment@conv) {

    if(integrate) {
      if(!all(is.na(est_early_rec_dev))) SE_Early <- sqrt(SD$diag.cov.random[names(SD$par.random) == "log_early_rec_dev"])
      SE_Main <- sqrt(SD$diag.cov.random[names(SD$par.random) == "log_rec_dev"])
    } else {
      if(!all(is.na(est_early_rec_dev))) SE_Early <- sqrt(diag(SD$cov.fixed)[names(SD$par.fixed) == "log_early_rec_dev"])
      SE_Main <- sqrt(diag(SD$cov.fixed)[names(SD$par.fixed) == "log_rec_dev"])
    }

    SE_Early2 <- est_early_rec_dev
    if(!all(is.na(est_early_rec_dev))) {
      SE_Early2[!is.na(SE_Early2)] <- SE_Early
    }
    SE_Main2 <- est_rec_dev
    SE_Main2[!is.na(SE_Main2)] <- SE_Main

    SE_Dev <- structure(c(rev(SE_Early2), SE_Main2), names = YearDev)

    first_non_zero <- which(!is.na(SE_Dev))[1]
    if(!is.na(first_non_zero) && first_non_zero > 1) {
      Dev_out <- Dev_out[-c(1:(first_non_zero - 1))]
      SE_Dev <- SE_Dev[-c(1:(first_non_zero - 1))]
      SE_Dev[is.na(SE_Dev)] <- 0
    }
    
    Assessment@Dev <- Dev_out
    Assessment@SE_Dev <- SE_Dev
  }
  
  
  ref_pt <- SCA_GTG_MSY_calc(report, info$data)
  report <- c(report, ref_pt)
  Assessment@FMSY <- report$FMSY
  Assessment@MSY <- report$MSY
  Assessment@BMSY <- report$BMSY
  Assessment@SSBMSY <- report$EMSY
  Assessment@VBMSY <- report$VBMSY
  Assessment@F_FMSY <- structure(report$F/report$FMSY, names = Year)
  Assessment@B_BMSY <- structure(report$B/report$BMSY, names = Yearplusone)
  Assessment@SSB_SSBMSY <- structure(report$E/report$EMSY, names = Yearplusone)
  Assessment@VB_VBMSY <- structure(report$VB/report$VBMSY, names = Yearplusone)
  Assessment@TMB_report <- report
  return(Assessment)
}
class(SCA_GTG) <- "Assess"



SCA_GTG_MSY_calc <- function(report, dat) {
  TMB_data <- list(model = "LeesApprox_MSY",
                   max_age = dat$max_age, M = dat$M, WAA = dat$WAA, mat = dat$mat,
                   SR_type = dat$SR_type, ngtg = dat$ngtg, Nbins = dat$Nbins,
                   LenBins = dat$LenBins, SAA = report$SAA, Select_at_length = report$Select_at_length,
                   LAA = dat$LAA, xout = dat$xout, distGTG = dat$distGTG, rdist = dat$rdist,
                   interp_check = dat$interp_check, interp_check2 = dat$interp_check2, integ_check = dat$integ_check,
                   integ_fac = dat$integ_fac, integ_ind = dat$integ_ind, use_LeesEffect = dat$use_LeesEffect)

  TMB_params <- list(log_F = log(0.1), Arec = report$Arec, Brec = report$Brec)
  map <- list()
  map$Arec <- map$Brec <- factor(NA)

  obj <- MakeADFun(TMB_data, TMB_params, map = map, DLL = "LeesApproxTMB", silent = TRUE)
  opt <- optimize(obj$fn, log(c(1e-8, 3)))
  report <- obj$report(obj$env$last.par.best)

  FMSY <- report$F
  MSY <- report$Yield
  VBMSY <- report$VB
  RMSY <- report$R
  BMSY <- report$B
  EMSY <- report$E

  return(list(FMSY = FMSY, MSY = MSY, VBMSY = VBMSY, RMSY = RMSY, BMSY = BMSY, EMSY = EMSY))
}

