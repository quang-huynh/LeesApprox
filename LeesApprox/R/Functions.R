addLines <- function(LenBins) {
  for (x in 1:length(LenBins)) {
    abline(v=LenBins[x], col="lightgray", lty=3)
  }
}

addBars <- function(probdf, col2) {
  yrs <- unique(probdf$Yr)
  temp <- probdf %>% dplyr::filter(Yr==yrs[2])
  for (x in 1:nrow(probdf)) {
    polygon(x=c(temp$Bin1[x], temp$Bin2[x], temp$Bin2[x], temp$Bin1[x]),
            y=c(0, 0, temp$Prob[x], temp$Prob[x]), col=col2)
  }
}

#' Title
#'
#' @param yr.st
#' @param yr.end
#' @param curF
#' @param F.pat
#' @param Fcv
#' @param plot
#'
#' @return
#' @export
#'
#' @examples
Ftrend <- function(yr.st, yr.end, curF,
                   F.pat= c('stable', 'inc', 'dec'),
                   Fcv=0.3, plot=TRUE) {
  F.pat <- match.arg(F.pat)
  yrs <- yr.st:yr.end
  yr.mid <- ceiling(mean(yrs))
  yr.ind <- which(yrs == yr.mid)
  nyrs <- length(yrs)
  Ferr <- exp(rnorm(nyrs, 0, Fcv))
  if (F.pat == "inc") {
    Ftrend <- seq(from=0, to=curF, length.out=nyrs) * Ferr
    Ftrend <- Ftrend/(Ftrend[length(Ftrend)]/curF)

  } else if (F.pat == "dec") {
    Ftrend <- c(seq(from=0, to=2*curF, length.out=yr.ind),
                seq(from=2*curF, to=curF, length.out=nyrs-yr.ind)) * Ferr
    Ftrend <- Ftrend/(Ftrend[length(Ftrend)]/curF)
  } else if (F.pat == "stable") {
    Ftrend <- c(seq(from=0, to=curF, length.out=yr.ind),
                seq(from=curF, to=curF, length.out=nyrs-yr.ind)) * Ferr
  }

  if (plot) {
    plot(yrs, Ftrend, type="l", bty="l", las=1, xlab="Years", ylab='apical Fishing mortality')
  }

  Ftrend
}

BHSRR <- function(SBcurr, SB0, R0, steepness) {
  (4 * R0 * steepness * SBcurr)/(SB0/R0 * R0 * (1-steepness) + (5*steepness-1)*SBcurr)
}

#' Title
#'
#' @param Linf
#' @param K
#' @param t0
#' @param M
#' @param L50
#' @param L95
#' @param LFS
#' @param L5
#' @param Vmaxlen
#' @param sigmaR
#' @param steepness
#' @param annualF
#' @param alpha
#' @param beta
#' @param LinfCV
#' @param ngtg
#' @param maxsd
#' @param binwidth
#'
#' @return
#' @export
#'
#' @examples
GTGpopsim <- function(Linf, K, t0, M, L50, L95, LFS, L5, Vmaxlen, sigmaR, steepness,
                      annualF,alpha=1E-5, beta=3, LinfCV=0.1, ngtg=101, maxsd=2, binwidth=1) {

  # Growth-type-groups
  distGTG <- seq(from=-maxsd, to=maxsd, length.out = ngtg)
  rdist <- dnorm(distGTG, 0, 1)/sum(dnorm(distGTG, 0, 1))
  Linfgtg <- Linf + LinfCV*Linf*distGTG

  LenBins <- seq(0, to=max(Linfgtg), by=binwidth)
  Nbins <- length(LenBins) -1
  By <- LenBins[2] - LenBins[1]
  LenMids <- seq(from=By*0.5, by=By, length.out = Nbins)

  nyrs <- length(annualF)

  maxage <- ceiling(-log(0.01)/M) # maxium age
  ages <- 1:maxage # in years

  maxAgeind <- maxage # maxage + 1
  ageVec <- 1:maxAgeind
  ind <- as.matrix(expand.grid(1:nyrs, ageVec,1:ngtg))
  LAA <- array(NA, dim=c(nyrs, maxAgeind, ngtg)) # Length-at-age by GTG and year
  LAA[ind] <- (Linfgtg[ind[,3]] * (1-exp(-K *(ages[ind[,2]]-t0))))
  WAA <- alpha * LAA^beta  # Weight-at-age by GTG and year

  # Maturity
  L50gtg <- L95gtg <- array(NA, dim=c(nyrs, maxAgeind, ngtg))
  L50gtg[ind] <- L50/Linf * Linfgtg  # assume constant L50/Linf but allow L50 to vary by year
  L95gtg[ind] <- L50gtg[ind] + (L95-L50)
  MAA <- 1/(1 + exp(-log(19) * ((LAA - L50gtg)/(L95gtg-L50gtg)))) # Maturity-at-age by GTG and year

  # selectivity-at-length - fishery
  sl <- (LFS - L5) /((-log(0.05,2))^0.5)
  sr <- (Linf - LFS) / ((-log(Vmaxlen,2))^0.5) # selectivity parameters are constant for all years
  SAA <- array(dnormal(LAA, LFS, sl, sr), dim=c(nyrs, maxAgeind, ngtg)) # selectivty-at-age by GTG and year
  SL <- dnormal(LenMids, LFS, sl, sr)

  M_array <- array(M, dim=c(nyrs, maxAgeind, ngtg))
  FAA <- array(NA, dim=c(nyrs, maxAgeind, ngtg))
  FAA[ind] <- SAA * annualF[ind[,1]] # fishing mortality at age by GTG and year
  ZAA <- FAA + M_array # Z-at-age by GTG and year

  R0 <- 1000
  # Unfished Year One
  Nunfished <- array(NA, dim=c(nyrs, maxAgeind, ngtg))
  SB <- array(NA, dim=c(nyrs, maxage, ngtg))
  Nunfished[1,1,] <- rdist * R0 # distribute virgin recruitment
  Nunfished[1,2:maxAgeind,] <- matrix(Nunfished[1,1,], nrow=maxage-1, ncol=ngtg, byrow=TRUE) *
    exp(-apply(M_array[1,ageVec-1,], 2, cumsum))

  SB[1,,] <- Nunfished[1,,] * WAA[1,,] * MAA[1,,]
  SB0 <- sum(SB[1,,])
  SBcurr <- Rec <- rep(NA, nyrs+1)
  SBcurr[1] <- SB0
  Rec[1] <- R0

  Nfished <- Nunfished
  recmu <- -0.5 * (sigmaR)^2
  recdevs <- exp(rnorm(nyrs, recmu, sigmaR))
  for (yr in 2:nyrs) {
    # message("Year ", yr, " of ", nyrs)
    Rec[yr] <- R0 # BHSRR(SBcurr[yr-1], SB0, R0, steepness) # recruitment
    Nfished[yr,1,] <- Rec[yr] * recdevs[yr] * rdist
    Nfished[yr,2:maxAgeind,] <- Nfished[yr-1,1:(maxage-1),] * exp(-ZAA[yr-1,1:(maxage-1),])
    SB[yr,,] <- Nfished[yr,,] * WAA[yr,,] * MAA[yr,,]
    SBcurr[yr] <- sum(SB[yr,,])
  }

  df <- data.frame(Yr=1:nyrs, Age=rep(ages, each=nyrs), N=as.vector(Nfished),
                   Select=as.vector(SAA), Length=as.vector(LAA),
                   GTG=rep(1:ngtg, each=(nyrs*(maxage))),
                   Linf=rep(Linfgtg, each=(nyrs*(maxage))))
  out <- list(df=df, LenBins=LenBins, LenMids=LenMids, N=Nfished)
  out
}


dnormal<-function(lens,lfs,sl,sr){
  cond<-lens<=lfs
  sel<-rep(NA,length(lens))
  sel[cond]<-2.0^-((lens[cond]-lfs)/sl*(lens[cond]-lfs)/sl)
  sel[!cond]<-2.0^-((lens[!cond]-lfs)/sr*(lens[!cond]-lfs)/sr)
  sel
}
