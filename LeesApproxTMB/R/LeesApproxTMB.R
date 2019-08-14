
#Sim <- LeesApprox(FVec=annualF, ngtg=ngtgVec[x], maxsd, binwidth=bw, M,
#                  Linf, K, t0,  LFS, L5, Vmaxlen, LinfCV, maxage)

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

  TMB_data <- list(FVec = FVec, ngtg = ngtg, LenBins = LenBins, LenMids = LenMids, ages = ages,
                   M = M, Linf = Linf, LAA = LAA, xout = xout, LFS = LFS, L5 = L5, Vmaxlen = Vmaxlen,
                   maxage = maxage, distGTG = distGTG, rdist = rdist)

  obj <- MakeADFun(data = TMB_data, parameters = list(p = 0), silent = TRUE, DLL = "LeesApproxTMB")
  report <- obj$report()

  return(list(probLA = report$probLA, LenComp = report$LenComp, Select_at_age = report$Select_at_age,
              Select_at_length = report$Select_at_length, LAA = report$LAA, LenMids = report$LenMids,
              LenBins = report$LenBins, NAA = report$NAA))
}

