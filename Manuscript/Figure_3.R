Figure3 <- function(pars, height=4, width=6, res=400, save=TRUE) {
  
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
  
  xx <- 1
  storeLAA <- storeSa <- NAA_store <- storepLAA <- list()
  DFlist <- list()
  for (Vmax in c(pars$Vmax, 1)) {
    for (curF in FVec) {
      F.trend <- rep(curF, length(yr.st: yr.end)) # Ftrend(yr.st, yr.end, curF, 'stable', plot=FALSE)
      
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
      storepLAA[[xx]] <- pLAA
      for (a in 1:maxage) Sa[a] <- sum(pLAA[a,] * SAL)
      storeSa[[xx]] <- Sa
      
      ## Length-at-age
      storeLAA[[xx]] <-  sapply(1:maxage, function(i) weighted.mean(LAA[i,], NAA[i,]))
      NAA_store[[xx]] <- list(LAA, NAA)
      DFlist[[xx]] <- data.frame(Ages=rep(1:24, 7), N=as.vector(NAA), L=as.vector(LAA),
                                 GTG=rep(1:7, each=24), Group=xx)
      xx <- xx + 1
    }
  }
  
  if (save) png('Figures/Figure3.png', units='in', height=height, width=width,
                res=res)
  
  par(mfrow=c(2,3), oma=c(2,2,1,1), mar=c(1,3,1,1))
  
  s_at_age <- do.call('cbind', storeSa)
  l_at_age <- do.call('cbind', storeLAA)
  w_at_age <- 1E-05*l_at_age^3
  
  Xat <- seq(0, maxage, by=4)
  xline <- 2
  yline <- 2.5
  textcex <- 0.7
  
  # Dome-shaped Selectivity
  sims <- 1:3
  matplot(s_at_age[,sims], type="l", bty="n", axes=FALSE, xlab="", ylab="",
          ylim=c(0, max(s_at_age)), yaxs="i", lwd=2, xpd=NA)
  axis(side=1, at=Xat, cex.axis=textcex, labels=FALSE)
  # mtext(side=1, line=xline, "Age", cex=textcex)
  axis(side=2, las=1, cex.axis=textcex)
  mtext(side=2, "Selectivity", cex=textcex, line=yline)
  text(1.5, max(s_at_age), "a)", xpd=NA, cex=textcex)
  
  legend('bottomright', legend=pars$FVec, col=1:3, lty=1:3, bty="n",
         title="Fishing mortality")
  
  matplot(l_at_age[,sims], type="l", bty="n", axes=FALSE, xlab="", ylab="",
          ylim=c(0, max(l_at_age)), yaxs="i", lwd=2, xpd=NA)
  axis(side=1, at=Xat, cex.axis=textcex, labels=FALSE)
  # mtext(side=1, line=xline, "Age", cex=textcex)
  axis(side=2, las=1, cex.axis=textcex)
  mtext(side=2, "Length", cex=textcex, line=yline)
  text(1.5, max(l_at_age), "b)", xpd=NA, cex=textcex)
  
  matplot(w_at_age[,sims], type="l", bty="n", axes=FALSE, xlab="", ylab="",
          ylim=c(0, max(w_at_age)), yaxs="i", lwd=2, xpd=NA)
  axis(side=1, at=Xat, cex.axis=textcex, labels=FALSE)
  # mtext(side=1, line=xline, "Age", cex=textcex)
  axis(side=2, las=1, cex.axis=textcex)
  mtext(side=2, "Weight", cex=textcex, line=yline)
  text(1.5, max(w_at_age), "c)", xpd=NA, cex=textcex)
  
  # Asymptotic Selectivity
  sims <- 4:6
  matplot(s_at_age[,sims], type="l", bty="n", axes=FALSE, xlab="", ylab="",
          ylim=c(0, max(s_at_age)), yaxs="i", lwd=2, xpd=NA)
  axis(side=1, at=Xat, cex.axis=textcex)
  mtext(side=1, line=xline, "Age", cex=textcex)
  axis(side=2, las=1, cex.axis=textcex)
  mtext(side=2, "Selectivity", cex=textcex, line=yline)
  text(1.5, max(s_at_age), "d)", xpd=NA, cex=textcex)

  matplot(l_at_age[,sims], type="l", bty="n", axes=FALSE, xlab="", ylab="",
          ylim=c(0, max(l_at_age)), yaxs="i", lwd=2, xpd=NA)
  axis(side=1, at=Xat, cex.axis=textcex)
  mtext(side=1, line=xline, "Age", cex=textcex)
  axis(side=2, las=1, cex.axis=textcex)
  mtext(side=2, "Length", cex=textcex, line=yline)
  text(1.5, max(l_at_age), "e)", xpd=NA, cex=textcex)
  
  matplot(w_at_age[,sims], type="l", bty="n", axes=FALSE, xlab="", ylab="",
          ylim=c(0, max(w_at_age)), yaxs="i", lwd=2, xpd=NA)
  axis(side=1, at=Xat, cex.axis=textcex)
  mtext(side=1, line=xline, "Age", cex=textcex)
  axis(side=2, las=1, cex.axis=textcex)
  mtext(side=2, "Weight", cex=textcex, line=yline)
  text(1.5,max(w_at_age), "f)", xpd=NA, cex=textcex)
  
  
  if (save) dev.off()
}
