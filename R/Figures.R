
#' Title
#'
#' @param height
#' @param width
#' @param res
#' @param Pars
#' @param save
#'
#' @return
#' @export
#'
#' @importFrom dplyr  %>%  filter group_by mutate summarize ungroup
Figure1 <- function(height=4, width=6, res=400, save=TRUE) {
  if (save) png('Figures/Figure1.png', units='in', height=height, width=width,
                res=res)

  binwidth <- 2
  ngtg_low <- 21
  ngtg_high <- 1001
  age <- 15
  yr <- c(1,68)
  par(mfcol=c(2,3), oma=c(4,4,0,0), mar=c(1,1,1,1))
  pch1 <- 21
  pch2 <- 24
  col1 <- 'black'
  col2 <- 'darkgray'
  xlab <- "Length class"
  ylab <- "Relative frequency"
  xline <- 2
  yline <- 2
  m.cex <- 1
  yaxs <- "i"
  xaxs <- "r"
  xpos <- function(df) min(df$Length)
  ypos <- function(maxy, ind=1) maxy$maxy[ind]
  tex.cex <- 0.8

    # Life-history parameters #
  M <- 0.2
  Linf <- 100
  LinfCV <- 0.1
  K <- 0.13333
  t0 <- 0
  sigmaR <- 0
  steepness <- 0.99

  L50 <- 66 # maturity
  L95 <- 70

  # Selectivity
  L5 <- 40
  LFS <- 50
  Vmaxlen <- 1

  maxsd <- 3

  # Trend in fishing mortality
  curF <- 0.4
  yr.st <- 1950
  yr.end <- 2017
  F.trend <- Ftrend(yr.st, yr.end, curF, 'stable', plot=FALSE)
  nyrs <- length(F.trend)

  # Low Res GTG #
  SimPop <- GTGpopsim(Linf, K, t0, M, L50, L95, LFS, L5, Vmaxlen, sigmaR,
                      steepness, F.trend, ngtg=ngtg_low, binwidth = binwidth,
                      maxsd=maxsd)

  LenMids <- SimPop$LenMids
  LenBins <- SimPop$LenBins
  t1 <- SimPop$df %>% filter(Yr %in% yr, Age==age)
  LenComp <- t1 %>% group_by(Yr, Age) %>% mutate(Bin=cut(Length, LenBins)) %>%
    group_by(Yr, Age, Bin) %>% summarize(P=sum(N)) %>% ungroup() %>%
    group_by(Yr, Age) %>% mutate(Prob = P/sum(P), Bin = as.character(Bin))

  n1 <- gsub("\\(", "", x=LenComp$Bin) %>% strsplit(",") %>%
    data.frame(stringsAsFactors=FALSE)
  Bin1 <- n1[1,] %>% as.numeric()
  Bin2 <- gsub("\\]", "", x=n1[2,]) %>% as.numeric()

  ProbDF <- data.frame(Yr=LenComp$Yr, Prob=LenComp$Prob, Bin1=Bin1, Bin2=Bin2)
  maxy <- ProbDF %>% group_by(Yr) %>% summarize(maxy=max(Prob))

  unfished <- t1 %>% filter(Yr==yr[1])
  fished <- t1 %>% filter(Yr==yr[2])
  plot(range(unfished$Length), c(0,1), ylim=c(0, max(ProbDF$Prob)), type="n",
       bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
  axis(side=1, labels=FALSE)
  axis(side=2)
  addLines(LenBins)
  points(unfished$Length, unfished$N/max(unfished$N)*maxy$maxy[1], pch=pch1,
         col="black", bg=col1, xpd=NA)
  points(fished$Length, fished$N/max(fished$N)*maxy$maxy[2], pch=pch2,
         col='black', bg=col2, xpd=NA)
  text(xpos(fished), ypos(maxy), "a)", cex=tex.cex, xpd=NA)

  plot(range(unfished$Length), c(0,1), ylim=c(0, maxy$maxy[2]), type="n",
       bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
  axis(side=1)
  axis(side=2)
  addLines(LenBins)
  addBars(ProbDF, col2)
  text(xpos(fished), ypos(maxy,2), "b)", cex=tex.cex, xpd=NA)
  mtext(side=2, line=yline, ylab, outer=TRUE, cex=m.cex, xpd=NA)

  # High Res GTG #
  SimPop <- GTGpopsim(Linf, K, t0, M, L50, L95, LFS, L5, Vmaxlen, sigmaR,
                      steepness, F.trend, ngtg=ngtg_high, binwidth = binwidth,
                      maxsd=maxsd)

  LenMids <- SimPop$LenMids
  LenBins <- SimPop$LenBins
  t1 <- SimPop$df %>% filter(Yr %in% yr, Age==age)
  LenComp <- t1 %>% group_by(Yr, Age) %>% mutate(Bin=cut(Length, LenBins)) %>%
    group_by(Yr, Age, Bin) %>% summarize(P=sum(N)) %>% ungroup() %>%
    group_by(Yr, Age) %>% mutate(Prob = P/sum(P), Bin = as.character(Bin))

  n1 <- gsub("\\(", "", x=LenComp$Bin) %>% strsplit(",") %>%
    data.frame(stringsAsFactors=FALSE)
  Bin1 <- n1[1,] %>% as.numeric()
  Bin2 <- gsub("\\]", "", x=n1[2,]) %>% as.numeric()

  ProbDF <- data.frame(Yr=LenComp$Yr, Prob=LenComp$Prob, Bin1=Bin1, Bin2=Bin2)
  maxy <- ProbDF %>% group_by(Yr) %>% summarize(maxy=max(Prob))

  unfished <- t1 %>% filter(Yr==yr[1])
  fished <- t1 %>% filter(Yr==yr[2])
  plot(range(unfished$Length), c(0,1), ylim=c(0, max(ProbDF$Prob)), type="n",
       bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
  axis(side=1, labels=FALSE)
  axis(side=2, labels=FALSE)
  addLines(LenBins)
  points(unfished$Length, unfished$N/max(unfished$N)*maxy$maxy[1], pch=pch1,
         col="black", bg=col1, xpd=NA)
  points(fished$Length, fished$N/max(fished$N)*maxy$maxy[2], pch=pch2,
         col='black', bg=col2, xpd=NA)
  mtext(side=1, line=xline, xlab, outer=TRUE, cex=m.cex, xpd=NA)
  text(xpos(fished), ypos(maxy), "c)", cex=tex.cex, xpd=NA)


  plot(range(unfished$Length), c(0,1), ylim=c(0, maxy$maxy[2]), type="n",
       bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
  axis(side=1)
  axis(side=2, labels=FALSE)
  addLines(LenBins)
  addBars(ProbDF, col2)
  text(xpos(fished), ypos(maxy,2), "d)", cex=tex.cex, xpd=NA)


  # Approx Method Res GTG
  SimPop <- GTGpopsim(Linf, K, t0, M, L50, L95, LFS, L5, Vmaxlen, sigmaR,
                      steepness, F.trend, ngtg=ngtg_low, binwidth = binwidth,
                      maxsd=maxsd)
  LenMids <- SimPop$LenMids
  LenBins <- SimPop$LenBins

  t1 <- SimPop$df %>% filter(Yr %in% yr, Age==age)
  xout <- sort(c(LenBins, t1$Length %>% unique()))

  # unfished & unfished
  Problist <- list()
  for (i in seq_along(yr)) {
    tempdf <- t1 %>% filter(Yr==yr[i])
    yout <- linear_int(tempdf$Length, tempdf$N, xout)
    Prob <- calcprob(x=tempdf$Length, y=tempdf$N, xout=xout, LenBins)
    Problist[[i]] <- data.frame(Yr=yr[i], Prob=Prob,
                                Bin1=LenBins[1:(length(LenBins)-1)],
                                Bin2=LenBins[2:(length(LenBins))])
  }
  ProbDF <- do.call("rbind", Problist)
  maxy <- ProbDF %>% group_by(Yr) %>% summarize(maxy=max(Prob))

  unfished <- t1 %>% filter(Yr==yr[1])
  fished <- t1 %>% filter(Yr==yr[2])
  plot(range(unfished$Length), c(0,1), ylim=c(0, max(ProbDF$Prob)), type="n",
       bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
  axis(side=1, labels=FALSE)
  axis(side=2, labels=FALSE)
  addLines(LenBins)
  lines(unfished$Length, unfished$N/max(unfished$N)*maxy$maxy[1], pch=pch1,
        col="black", bg=col1, type="b", xpd=NA)
  lines(fished$Length, fished$N/max(fished$N)*maxy$maxy[2], pch=pch2,
        col='black', bg=col2, type="b", xpd=NA)
  text(xpos(fished), ypos(maxy), "e)", cex=tex.cex, xpd=NA)

  plot(range(unfished$Length), c(0,1), ylim=c(0, maxy$maxy[2]), type="n",
       bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
  axis(side=1)
  axis(side=2, labels=FALSE)
  addLines(LenBins)
  addBars(ProbDF, col2)
  text(xpos(fished), ypos(maxy,2), "f)", cex=tex.cex, xpd=NA)

  if(save) dev.off()
}

#' Title
#'
#' @param height
#' @param width
#' @param res
#' @param save
#'
#' @return
#' @export
#'
#' @examples
Figure2 <- function(height=2.5, width=6, res=400, save=TRUE) {

  ngtg <- 7
  LinfCV <- 0.1
  Linf <- 100
  maxsd <- 2
  binwidth <- 5

  K <- 0.25
  t0 <- 0
  M <- 0.2
  maxage <- ceiling(-log(0.01)/M)
  Ages <- 1:maxage

  L5 <- 10
  LFS <- 30
  Vmax <- 0.6

  curF <- 0.2
  yr.st <- 1950
  yr.end <- 2017
  F.trend <- Ftrend(yr.st, yr.end, curF, 'stable', plot=FALSE)
  Curr.Yr <- length(F.trend)
  HistFs <- F.trend[(Curr.Yr - maxage+1): Curr.Yr]

  distGTG <- seq(-maxsd, maxsd, length.out=ngtg)
  Linfgtg <- Linf + LinfCV*Linf*distGTG
  LenBins <- seq(0, to=max(Linfgtg), by=binwidth)
  Nbins <- length(LenBins) -1
  By <- LenBins[2] - LenBins[1]
  LenMids <- seq(from=By*0.5, by=By, length.out = Nbins)

  VB <- function(Linf, K, t0, age) Linf * (1-exp(-K*(age-t0)))
  LAA <- lapply(1:ngtg, function(x) VB(Linfgtg[x], K, t0, Ages))
  LAA <- do.call("cbind", LAA)

  Sigma_asc <- (L5-LFS)/sqrt(-log(0.05,2))
  Sigma_desc <- (Linf-LFS)/sqrt(-log(Vmax,2))

  SAA <- lapply(1:ngtg, function(x) dnormal(LAA[,x],LFS, Sigma_asc, Sigma_desc))
  SAA <- do.call("cbind", SAA)
  SAL <- dnormal(LenMids, LFS, Sigma_asc, Sigma_desc)

  pg <- dnorm(distGTG, 0 ,1); pg <- pg/sum(pg)

  NAA <- matrix(NA, nrow=maxage, ncol=ngtg)
  NAA[1, ] <- pg

  for (a in 2:maxage) {
    Z <- apply(matrix(HistFs[(length(HistFs)-1) : (length(HistFs)-a+1)], nrow=a-1, ncol=ngtg) * SAA[1:(a-1), ], 2, sum) + M*(a-1)
    NAA[a,] <- pg * exp(-Z)
  }

  if (save) png('Figures/Figure2.png', units='in', height=height, width=width, res=res)
  par(mfrow=c(1,3), oma=c(2,2,0,0), mar=c(2,2,2,2))

  Xat <- seq(0, maxage, by=4)
  xline <- 2
  yline <- 2.5
  textcex <- 0.7
  matplot(LAA, type="b", bty="n", axes=FALSE, xlab="", ylab="",
          ylim=c(0, max(LAA)), yaxs="i", lwd=2)
  axis(side=1, at=Xat, cex.axis=textcex)
  mtext(side=1, line=xline, "Age", cex=textcex)
  axis(side=2, las=1, cex.axis=textcex)
  mtext(side=2, "Length", cex=textcex, line=yline)
  text(1, max(LAA), "a)", xpd=NA, cex=textcex)

  matplot(SAA, type="b", bty="n", axes=FALSE, xlab="", ylab="",
          ylim=c(0, 1.025), yaxs="i", lwd=2, xpd=NA)
  axis(side=1, at=Xat, cex.axis=textcex)
  mtext(side=1, line=xline, "Age", cex=textcex)
  axis(side=2, las=1, cex.axis=textcex)
  mtext(side=2, "Selectivity", cex=textcex, line=yline)
  text(1, 1.02, "b)", xpd=NA, cex=textcex)

  age <- 15

  xs <- LAA[age, ]
  ys <- NAA[age, ]/sum(NAA[age, ])
  xout <- sort(c(xs, LenBins))
  Yout <- linear_int(xs, ys, xout)

  cond <- xout %in% LenBins
  ind <- which(Yout>0)
  ind <- c(ind[1]-1, ind, ind[length(ind)]+1)
  plot(xout[ind], Yout[ind], type="n", ylim=c(0, 0.3), yaxs="i", axes=FALSE,
       xlab="", ylab="", las=1)
  axis(side=1, cex.axis=textcex)
  mtext(side=1, line=2, "Length class", cex=textcex)
  axis(side=2, las=1, cex.axis=textcex)
  mtext(side=2, line=2.5, "Relative Frequency", cex=textcex)
  lines(xout[!cond], Yout[!cond], type="b", pch=15, xpd=NA)
  # text(xout[!cond], Yout[!cond]+.0125, 1:ngtg, pch=15, cex=1.5, xpd=NA)
  ind1 <- which(Yout[cond]>0)
  ind1 <- c(ind1[1]-1, ind1, ind1[length(ind1)]+1)
  points(xout[cond][ind1], Yout[cond][ind1], pch=16, xpd=NA)
  abline(v=LenBins, lty=2, col="lightgray")
  text(75, 0.3, "c)", xpd=NA, cex=textcex)

  ind <- which(xout==90)
  density <- 20
  polygon(x=c(90, 95, 95, 90),
          y=c(0, 0, Yout[ind],Yout[ind]), border=1, lty=1, density =density)

  polygon(x=c(90, xout[ind+1], xout[ind+1], 90),
          y=c(Yout[ind], Yout[ind], Yout[ind+1], Yout[ind]),
          border=1, lty=1, density=density, angle=-45)

  polygon(x=c(xout[ind+1], 95, 95, xout[ind+1]),
          y=c(Yout[ind], Yout[ind], Yout[ind+2], Yout[ind+1]),
          border=1, lty=1, density=density, angle=-45)
  polygon(x=c(xout[ind+1], 95, 95, xout[ind+1]),
          y=c(Yout[ind], Yout[ind], Yout[ind+2], Yout[ind+1]),
          border=1, lty=1, density=density, angle=45)

  if (save) dev.off()
}

#' Title
#'
#' @param height
#' @param width
#' @param res
#' @param save
#'
#' @return
#' @export
#'
#' @examples
Figure3 <- function(height=2.5, width=6, res=400, save=TRUE) {

  ngtg <- 7
  LinfCV <- 0.1
  Linf <- 100
  maxsd <- 2
  binwidth <- 5

  K <- 0.25
  t0 <- 0
  M <- 0.2
  maxage <- ceiling(-log(0.01)/M)
  Ages <- 1:maxage

  L5 <- 10
  LFS <- 30
  Vmax <- 0.6

  yr.st <- 1950
  yr.end <- 2017
  FVec <- c(0, 0.2, 0.6)
  xx <- 1
  storeLAA <- storeSa <- list()
  for (curF in FVec) {
    F.trend <- Ftrend(yr.st, yr.end, curF, 'stable', plot=FALSE)

    Curr.Yr <- length(F.trend)
    HistFs <- F.trend[(Curr.Yr - maxage+1): Curr.Yr]

    distGTG <- seq(-maxsd, maxsd, length.out=ngtg)
    Linfgtg <- Linf + LinfCV*Linf*distGTG

    LenBins <- seq(0, to=max(Linfgtg), by=binwidth)
    Nbins <- length(LenBins) -1
    By <- LenBins[2] - LenBins[1]
    LenMids <- seq(from=By*0.5, by=By, length.out = Nbins)

    VB <- function(Linf, K, t0, age) Linf * (1-exp(-K*(age-t0)))
    LAA <- lapply(1:ngtg, function(x) VB(Linfgtg[x], K, t0, Ages))
    LAA <- do.call("cbind", LAA)

    Sigma_asc <- (L5-LFS)/sqrt(-log(0.05,2))
    Sigma_desc <- (Linf-LFS)/sqrt(-log(Vmax,2))

    SAA <- lapply(1:ngtg, function(x) dnormal(LAA[,x],LFS, Sigma_asc, Sigma_desc))
    SAA <- do.call("cbind", SAA)

    SAL <- dnormal(LenMids, LFS, Sigma_asc, Sigma_desc)

    pg <- dnorm(distGTG, 0 ,1); pg <- pg/sum(pg)
    NAA <- matrix(NA, nrow=maxage, ncol=ngtg)
    NAA[1, ] <- pg

    for (a in 2:maxage) {
      Z <- apply(matrix(HistFs[(length(HistFs)-1) : (length(HistFs)-a+1)], nrow=a-1, ncol=ngtg) * SAA[1:(a-1), ], 2, sum) + M*(a-1)
      NAA[a,] <- pg * exp(-Z)
    }

    # Calc Prob Length-at-Age
    pLAA <- matrix(NA, nrow=maxage, ncol=length(LenMids))
    for (age in 1:maxage) {
      xs <- LAA[age, ]
      ys <- NAA[age, ]/sum(NAA[age, ])
      xout <- sort(c(xs, LenBins))
      pLAA[age,] <- calcprob(xs, ys, xout, LenBins)
    }
    ages <- seq(2, maxage, by=4)

    # Selectivity-at-age
    Sa <- rep(NA, maxage)
    for (a in 1:maxage) Sa[a] <- sum(pLAA[a,] * SAL)
    storeSa[[xx]] <- Sa

    ## Length-at-age
    storeLAA[[xx]] <- apply(LAA, 1, weighted.mean, w=apply(NAA, 2, sum))
    xx <- xx + 1
  }
  if (save) png('Figures/Figure3.png', units='in', height=height, width=width,
                res=res)

  par(mfrow=c(1,3), oma=c(2,2,1,1), mar=c(1,3,1,1))

  s_at_age <- do.call('cbind', storeSa)
  l_at_age <- do.call('cbind', storeLAA)
  w_at_age <- 1E-05*l_at_age^3

  Xat <- seq(0, maxage, by=4)
  xline <- 2
  yline <- 2.5
  textcex <- 0.7

  matplot(s_at_age, type="l", bty="n", axes=FALSE, xlab="", ylab="",
          ylim=c(0, 1.025), yaxs="i", lwd=2, xpd=NA)
  axis(side=1, at=Xat, cex.axis=textcex)
  mtext(side=1, line=xline, "Age", cex=textcex)
  axis(side=2, las=1, cex.axis=textcex)
  mtext(side=2, "Selectivity", cex=textcex, line=yline)
  text(1.5, 1.02, "a)", xpd=NA, cex=textcex)


  matplot(l_at_age, type="l", bty="n", axes=FALSE, xlab="", ylab="",
          ylim=c(0, Linf), yaxs="i", lwd=2, xpd=NA)
  axis(side=1, at=Xat, cex.axis=textcex)
  mtext(side=1, line=xline, "Age", cex=textcex)
  axis(side=2, las=1, cex.axis=textcex)
  mtext(side=2, "Length", cex=textcex, line=yline)
  text(1.5, Linf, "b)", xpd=NA, cex=textcex)

  matplot(w_at_age, type="l", bty="n", axes=FALSE, xlab="", ylab="",
          ylim=c(0, 1E-5*Linf^3), yaxs="i", lwd=2, xpd=NA)
  axis(side=1, at=Xat, cex.axis=textcex)
  mtext(side=1, line=xline, "Age", cex=textcex)
  axis(side=2, las=1, cex.axis=textcex)
  mtext(side=2, "Weight", cex=textcex, line=yline)
  text(1.5, 1E-5*Linf^3, "c)", xpd=NA, cex=textcex)

  if (save) dev.off()
}
