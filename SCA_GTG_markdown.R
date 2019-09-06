library(MSEtool)

summary_SCA_GTG <- function(Assessment) {
  assign_Assessment_slots(Assessment)

  if(conv) current_status <- c(F_FMSY[length(F_FMSY)], B_BMSY[length(B_BMSY)], B_B0[length(B_B0)])
  else current_status <- c(NA, NA, B_B0[length(B_B0)])
  current_status <- data.frame(Value = current_status)
  rownames(current_status) <- c("F/FMSY", "B/BMSY", "B/B0")

  Value <- c(h, info$data$M[1], info$data$max_age, info$LH$Linf, info$LH$K, info$LH$t0,
             info$LH$a * info$LH$Linf ^ info$LH$b, info$LH$A50, info$LH$A95,
             info$data$ngtg, info$LH$LenCV)
  Description = c("Stock-recruit steepness", "Natural mortality", "Maximum age (plus-group)", "Asymptotic length", "Growth coefficient",
                  "Age at length-zero", "Asymptotic weight", "Age of 50% maturity", "Age of 95% maturity",
                  "Number of growth-type-groups", "Variability in length-at-age")
  rownam <- c("h", "M", "maxage", "Linf", "K", "t0", "Winf", "A50", "A95", "GTG", "LenCV")
  input_parameters <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
  rownames(input_parameters) <- rownam
  if(!"transformed_h" %in% names(obj$env$map)) input_parameters <- input_parameters[-1, ]

  if(conv) Value <- c(SSB0, MSY, FMSY, VBMSY, SSBMSY, SSBMSY/SSB0)
  else Value <- rep(NA, 6)

  Description <- c("Unfished spawning stock biomass (SSB)", "Maximum sustainable yield (MSY)", "Fishing mortality at MSY",
                   "Vulnerable biomass at MSY", "SSB at MSY", "Spawning depletion at MSY")
  derived <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
  rownames(derived) <- c("SSB0", "MSY", "FMSY", "VBMSY", "SSBMSY", "SSBMSY/SSB0")

  if(!is.character(SD)) {
    model_estimates <- summary(SD)[rownames(summary(SD)) != "log_rec_dev" & rownames(summary(SD)) != "log_early_rec_dev" &
                                     rownames(summary(SD)) != "logF", ]
    model_estimates <- model_estimates[!is.na(model_estimates[, 2]) && model_estimates[, 2] > 0, ]
    if(length(SE_Dev) == 0) SE_Dev <- rep(NA, length(Dev))
    dev_estimates <- cbind(Dev, SE_Dev)
    rownames(dev_estimates) <- paste0("log_rec_dev_", names(Dev))

    model_estimates <- rbind(model_estimates, dev_estimates)
  } else {
    model_estimates <- SD
  }

  output <- list(model = "Statistical Catch-at-Age with Growth Type Groups (SCA_GTG)",
                 current_status = current_status, input_parameters = input_parameters,
                 derived_quantities = derived, model_estimates = model_estimates)
  return(output)
}

rmd_SCA_GTG <- function(Assessment) {
  filename <- file.path(system.file(package = "LeesApproxTMB"), "report_SCA_GTG.Rmd")
  report <- readLines(filename)
  return(report)
}

plot_yield_SCA_GTG <- function(dat, report, fmsy, msy, xaxis = c("F", "Biomass", "Depletion")) {
  xaxis <- match.arg(xaxis)
  F.vector = seq(0, 2.5 * fmsy, length.out = 1e2)
  F.vector[1] <- 1e-8

  BMSY <- report$EMSY
  B0 <- report$E0

  TMB_data <- list(model = "LeesApprox_MSY",
                   max_age = dat$max_age, M = dat$M, WAA = dat$WAA, mat = dat$mat,
                   SR_type = dat$SR_type, ngtg = dat$ngtg, Nbins = dat$Nbins,
                   LenBins = dat$LenBins, SAA = report$SAA, Select_at_length = report$Select_at_length,
                   LAA = dat$LAA, xout = dat$xout, distGTG = dat$distGTG, rdist = dat$rdist,
                   use_LeesEffect = dat$use_LeesEffect)

  TMB_params <- list(log_F = log(0.1), Arec = report$Arec, Brec = report$Brec)
  map <- list()
  map$Arec <- map$Brec <- factor(NA)

  obj <- MakeADFun(TMB_data, TMB_params, map = map, DLL = "LeesApproxTMB", silent = TRUE)

  report <- lapply(log(F.vector), obj$report)

  Yield <- vapply(report, getElement, numeric(1), "Yield")
  R <- vapply(report, getElement, numeric(1), "R")
  Biomass <- vapply(report, getElement, numeric(1), "E")

  ind <- R >= 0

  if(xaxis == "F") {
    plot(F.vector[ind], Yield[ind], typ = 'l', xlab = "Fishing Mortality",
         ylab = "Equilibrium yield")
    segments(x0 = fmsy, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = fmsy, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if(xaxis == "Biomass") {
    plot(Biomass[ind], Yield[ind], typ = 'l', xlab = "Spawning Stock Biomass",
         ylab = "Equilibrium yield")
    segments(x0 = BMSY, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = BMSY, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if(xaxis == "Depletion") {
    plot(Biomass[ind]/B0, Yield[ind], typ = 'l',
         xlab = expression(SSB/SSB[0]), ylab = "Equilibrium yield")
    segments(x0 = BMSY/B0, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = BMSY/B0, lty = 2)
    abline(h = 0, col = 'grey')
  }
  invisible(data.frame(F = F.vector[ind], Yield = Yield[ind], B = Biomass[ind], B_B0 = Biomass[ind]/B0))
}


environment(summary_SCA_GTG) <- environment(rmd_SCA_GTG) <- environment(plot_yield_SCA_GTG) <-
  asNamespace("MSEtool")
