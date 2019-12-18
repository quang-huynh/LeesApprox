dnormal<-function(lens,lfs,sl,sr){
  cond<-lens<=lfs
  sel<-rep(NA,length(lens))
  sel[cond]<-2.0^-((lens[cond]-lfs)/sl*(lens[cond]-lfs)/sl)
  sel[!cond]<-2.0^-((lens[!cond]-lfs)/sr*(lens[!cond]-lfs)/sr)
  sel
}

calcprobR <- function(tempLens, tempNs, xout, LenBins) {
  # yout <- linear_int(tempLens, tempNs, xout);
  yout <- approx(tempLens, tempNs, xout)[[2]]
  yout[is.na(yout)] <- 0

  nclasses <- length(LenBins)-1
  Prob <- rep(0, nclasses)
  for (i in 1:nclasses) {
    L1 <- LenBins[i]
    L2 <- LenBins[i+1]
    ind <- (xout >=L1) & (xout <= L2)
    if (sum(ind)==1) {
      Prob[i] <- youtp[i]
    } else {
      Ys <- yout[ind]
      Xs <- xout[ind]
      ns <- sum(ind)-1
      area <- rep(0, ns)
      for (j in 1:ns) {
        p1 <- Ys[j]
        p2 <- Ys[j+1]
        h1 <- abs(p1-p2)
        by <- Xs[j+1] - Xs[j]
        tvec <- c(p1,p2)
        h2 <- min(tvec)
        area[j] <- h2 * by + 0.5 * by * h1
      }
      Prob[i] <- sum(area)
    }
  }
  Probout <- Prob/sum(Prob)
  Probout
}

LeesApproxR <- function(FVec, ngtg=1001, maxsd, binwidth, M,
                        Linf, K, t0,  LFS, L5, Vmaxlen, LinfCV, maxage) {

  distGTG <- seq(-maxsd, maxsd, length.out=ngtg)
  rdist <- dnorm(distGTG, 0, 1.0, 0)
  rdist <- rdist/sum(rdist)
  Linfgtg <- Linf + LinfCV*Linf*distGTG

  LenBins <- seq(0, max(Linfgtg), binwidth)
  LenMids <- LenBins - 0.5 * binwidth
  LenMids <- LenMids[LenMids>0]

  LAA <- matrix(NA, maxage, ngtg)
  ages <- 1:maxage

  for (g in 1:ngtg) {
    LAA[,g] <- Linfgtg[g] * (1-exp(-K*(ages-t0)))
  }

  # calculate selectivity-at-age for each GTG
  sl <- (LFS - L5) / sqrt(-log2(0.05))
  sr <- (Linf - LFS) / sqrt(-log2(Vmaxlen))

  SAA <- matrix(NA, maxage, ngtg)
  for (g in 1:ngtg) {
    SAA[,g] <-  dnormal(LAA[,g], LFS, sl, sr)
  }

  # calculate selectivity-at-length
  Select_at_length <- dnormal(LenMids, LFS, sl, sr);

  # create vector of Fs for last maxage years
  # F assumed 0 for years where no F are provided
  FVec2 <- rep(0, maxage)
  ll <- length(FVec)
  FVec2 <- FVec[(ll-maxage+1):ll]


  # loop over ages and calculate N per recruit for each age class
  Ns <- matrix(NA, maxage, ngtg)
  Ns[1,] <- rdist

  for (age in 2:maxage) {
    yr_st <- maxage-age
    age_end <- age-1
    Zs <- rep(0, ngtg)
    for (age2 in 1:age_end) {
      Zs <- Zs + M + FVec2[yr_st+age2] *SAA[age2,]
    }
    Ns[age,] <- Ns[1,] * exp(-Zs)
  }

  # Calculate prob L|A
  Nbins <- length(LenMids)
  probLA <- matrix(0, maxage, Nbins)
  for (age in 1:maxage) {
    tempLens <- LAA[age,]
    tempNs <- Ns[age,]
    xout <- c(LenBins, tempLens) %>% sort()
    probLA[age,] <- calcprobR(tempLens, tempNs, xout, LenBins)
  }

  # Calculate mean selectivity-at-age
  Select_at_age <- rep(0,maxage)
  for (age in 1:maxage) {
    Select_at_age[age] <- Select_at_age[age] + sum(probLA[age,] * Select_at_length)
  }

  # Calculate mean length-at-age
  Len_at_age <- rep(0, maxage)
  for (age in 1:maxage) {
    Len_at_age[age] <- weighted.mean(LAA[age,], Ns[age,]) # sum(LAA[age,] * Ns[age,])/sum(Ns[age,])
  }

  # Calculate fished length composition
  NAA <- apply(Ns, 1, sum)

  LenComp <- rep(0, Nbins)
  for (l in 1:Nbins) {
    LenComp[l] <- sum(probLA[,l] * Select_at_length[l] * NAA/sum(NAA))
  }
  LenComp = LenComp/sum(LenComp);

  out <- list()
  out[[1]] = probLA;
  out[[2]] = LenComp;
  out[[3]] = Select_at_age;
  out[[4]] = Select_at_length;
  out[[5]] = LAA;
  out[[6]] = LenMids;
  out[[7]] = LenBins;
  out[[8]] = NAA;
  out
}



