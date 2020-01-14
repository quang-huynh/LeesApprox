GenerateData <- function(Stock=1, DatYears=10, CobCV=0.1, IobCV=0.1,
                         LH_CV=0.05,
                         LengthSampSize=250,
                         AgeSampSize=250,
                         Fmulti=1, nbins=30) {
  
  st <- Stock
  StockNames <- c('Queen triggerfish',
                  'Stoplight parrotfish',
                  'Yellowtail snapper')
  
  Stocks <- list(Queen_Triggerfish_STT_NOAA,
                 Stoplight_Parrotfish_STX_NOAA,
                 Yellowtail_Snapper_PR_NOAA)
  message(StockNames[st])
  
  Stock <- Stocks[[st]]
  Linf <- mean(Stock@Linf)
  K <- mean(Stock@K)
  t0 <- mean(Stock@t0)
  M <- mean(Stock@M)
  maxage <- ceiling(-log(0.001)/M)
  
  L50 <- mean(Stock@L50)
  L95 <- L50 +  mean(Stock@L50_95)
  
  L5 <- mean(Stock@L5) * L50
  LFS <- mean(Stock@LFS) * L50
  Vmaxlen <- mean(Stock@Vmaxlen)
  sigmaR <- mean(Stock@Perr)
  steepness <- mean(Stock@h)
  
  alpha <- Stock@a
  beta <- Stock@b
  LinfCV <- 0.1
  
 
  
  LHpars <- data.frame(Stock=StockNames[st], Linf=Linf, K=K,
                       t0=t0, M=M, maxage=maxage, L50=L50, L95=L95,
                       L5=L5, LFS=LFS, Vmaxlen=Vmaxlen, sigmaR=sigmaR,
                       steepness=steepness, alpha=alpha, beta=beta, LinfCV=LinfCV)
  
  maxsd <- 2
  
  maxL <- Linf + LinfCV*Linf*maxsd
  
  bins <- seq(0, to=max(maxL), length.out=nbins+1)
  
  binwidth <- round((bins[2] - bins[1])/5) * 5
  
  Fm <- Fmulti * M
  annualF <- Ftrend(1970, 2019, Fm, 'stable', Fcv=0.1, plot=FALSE)
  annualF[1] <- annualF[2] * 0.5
  
  SimPop <- GTGpopsim(Linf, K, t0, M, L50, L95, LFS, L5, Vmaxlen, sigmaR,
                      steepness, annualF,alpha, beta, LinfCV, ngtg=1001, maxsd,
                      binwidth)
  
  DF <- SimPop[[1]]
  DF <- DF %>% group_by(Yr) %>% mutate(Catch=sum(CAA*Weight),
                                       TotalB=sum(N*Weight))
  
  TSData <- DF %>% select(Yr, Catch, Index=TotalB) %>% distinct()
  TSData$Index <- TSData$Index/mean(TSData$Index)
  
  # Catch-at-age
  DF$vulnN <- DF$N * DF$Select
  CAA_DF <- DF %>% group_by(Yr, Age) %>% summarize(CAA=sum(CAA))
  VAge_Comp <- CAA_DF %>% tidyr::pivot_wider(names_from='Yr',
                                             values_from="CAA")
  
  AgeSamps <- sapply(1:length(annualF), function(i) 
    rmultinom(n=1, size=AgeSampSize, VAge_Comp[[i+1]]))

  CAA <- t(AgeSamps)
  
  # Catch-at-length
  CAL <- matrix(0, nrow=length(annualF), ncol=length(SimPop$LenMids))
  for (yr in 1:length(annualF)) {
    df <- DF %>% filter(Yr==yr)
    lenP <- rep(0, length(SimPop$LenMids))
    for (l in seq_along(SimPop$LenMids)) {
      ind <- df$Length >= SimPop$LenBins[l] & df$Length < SimPop$LenBins[l+1]
      lenP[l] <- sum(df$vulnN[ind])
    }
    lenP <- lenP/sum(lenP)
    if (!all(lenP == 0)) {
      CAL[yr,] <- t(rmultinom(n=1, size=LengthSampSize, prob=lenP))
    } 
  }
  

  Cobs <- rlnorm(length(annualF), 0, CobCV)
  Iobs <- rlnorm(length(annualF), 0, IobCV)
  
  CatchDat <- TSData$Catch * Cobs
  IndexDat <- TSData$Index * Iobs
  
  if (!is.na(DatYears)) {
    AllYears <- 1:length(annualF)
    NAYears <- 1:(length(annualF) - DatYears)
    
    IndexDat[NAYears] <- NA
    CAL[NAYears,] <- NA
    CAA[NAYears,] <- NA 
  }

  
  Data <- new("Data")
  Data@CAL_bins <- SimPop$LenBins
  Data@CAL_mids <- SimPop$LenMids
  Data@CAL <- array(CAL, dim=c(1, length(annualF), length(SimPop$LenMids)))
  Data@CAA <- array(CAA, dim=c(1, length(annualF), max(DF$Age)))
  Data@Year <- unique(DF$Yr)
  Data@Cat <- matrix(CatchDat, nrow=1)
  Data@Ind <- matrix(IndexDat, nrow=1)
  
  Data@MaxAge <- max(DF$Age)
  Data@Mort <- M * rlnorm(1, 0, LH_CV)
  Data@vbt0 <- t0 * rlnorm(1, 0, LH_CV)
  Data@vbK <- K * rlnorm(1, 0, LH_CV)
  Data@vbLinf <- Linf * rlnorm(1, 0, LH_CV)
  Data@L50 <- L50 * rlnorm(1, 0, LH_CV)
  L50_95 <- L95 - L50 
  Data@L95 <- L50 + L50_95 * rlnorm(1, 0, LH_CV)
  Data@wla <- alpha
  Data@CV_vbLinf <- LinfCV
  Data@wlb <- beta
  Data@steep <- steepness
  Data@sigmaR <- sigmaR
  
  Data@LFS <- LFS
  Data@LFC <- L5
  
  Data@AddInd <- array(NA, dim=c(1,1,1))
  Data@CV_AddInd <- array(NA, dim=c(1,1,1))
  Data@AddIndV <- array(NA, dim=c(1,1,1))
  
  out <- list()
  out$Catch <- TSData$Catch
  out$Index <- TSData$Index
  out$annualF <- annualF 
  out$Data <- Data
  out$LHpars <- LHpars
  out
  
  
}