Figure2 <- function(pars,height=2.5, width=6, res=400, save=TRUE) {
  
  ngtg  <- pars$ngtg 
  LinfCV <- pars$LinfCV
  Linf <- pars$Linf
  maxsd  <- pars$maxsd 
  binwidth  <- pars$binwidth 
  K  <- pars$K 
  t0 <- pars$t0
  M  <- pars$M 
  maxage  <- pars$maxage 
  Ages  <- pars$Ages 
  L5  <- pars$L5 
  LFS  <- pars$LFS 
  Vmax <- pars$Vmax
  yr.st  <- pars$yr.st 
  yr.end  <- pars$yr.end 
  FVec <- pars$FVec 
  
  # F.trend <- Ftrend(yr.st, yr.end, curF, 'stable', plot=FALSE)
  F.trend <- rep(0, length(yr.st:yr.end))
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
    Z <- apply(matrix(HistFs[(length(HistFs)-1):(length(HistFs)-a+1)], 
                      nrow=a-1, ncol=ngtg) * SAA[1:(a-1), ], 2, sum) + M*(a-1)
    NAA[a,] <- pg * exp(-Z)
  }
  
  if (save) png('Figures/Figure2.png', units='in', height=height, width=width, res=res)
  par(mfrow=c(1,3), oma=c(2,2,0,0), mar=c(2,2,2,2))
  
  Xat <- seq(0, maxage, by=4)
  xline <- 2
  yline <- 2.5
  textcex <- 0.7
  matplot(LAA, type="b", bty="n", axes=FALSE, xlab="", ylab="",
          ylim=c(0, max(LAA)*1.02), yaxs="i", lwd=2, xpd=NA)
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
  
  minind <- which(xout == min(xout[!cond]))
  maxind <- which(xout == max(xout[!cond]))
  lines(c(xout[minind-1],xout[minind]), c(Yout[minind-1],Yout[minind]))
  lines(c(xout[maxind],xout[maxind+1]), c(Yout[maxind],Yout[maxind+1]))
  
  ind1 <- which(Yout[cond]>0)
  ind1 <- c(ind1[1]-1, ind1, ind1[length(ind1)]+1)
  points(xout[cond][ind1], Yout[cond][ind1], pch=16, xpd=NA)
  abline(v=LenBins, lty=2, col="lightgray")
  text(75, 0.3, "c)", xpd=NA, cex=textcex)
  
  text(xout[!cond], Yout[!cond]+.015, 1:ngtg, pch=15, cex=0.7, xpd=NA)
  
  ind <- which(xout==90)
  density <- 20
  border <- FALSE
  polygon(x=c(90, 95, 95, 90),
          y=c(0, 0, Yout[ind+1],Yout[ind]), border=border, lty=1, density =density)
  
  polygon(x=c(90, xout[ind+1], xout[ind+1], 90),
          y=c(Yout[ind], Yout[ind], Yout[ind+1], Yout[ind]),
          border=border, lty=1, density=density, angle=45)
  # 
  polygon(x=c(xout[ind+1], 95, 95, xout[ind+1]),
          y=c(Yout[ind], Yout[ind], Yout[ind+2], Yout[ind+1]),
          border=border, lty=1, density=density, angle=45)
  # polygon(x=c(xout[ind+1], 95, 95, xout[ind+1]),
  #         y=c(Yout[ind], Yout[ind], Yout[ind+2], Yout[ind+1]),
  #         border=1, lty=1, density=density, angle=45)
  
  if (save) dev.off()
}