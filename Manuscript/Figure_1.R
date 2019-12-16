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