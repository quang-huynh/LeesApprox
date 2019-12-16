StockNames <- c('Queen triggerfish',
                'Stoplight parrotfish',
                'Yellowtail snapper')

Stocks <- list(Queen_Triggerfish_STT_NOAA,
               Stoplight_Parrotfish_STX_NOAA,
               Yellowtail_Snapper_PR_NOAA)

xx <- 0 
DFstore <- list()

for (st in 1:length(Stocks)) {
  message(StockNames[st])
  Stock <- Stocks[[st]]
  Linf <- mean(Stock@Linf)
  K <- mean(Stock@K)
  t0 <- mean(Stock@t0)
  M <- mean(Stock@M)
  maxage <- ceiling(-log(0.01)/M)
  
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
  
  maxsd <- 2
  
  bwVec <- c(2,5,10)
  ngtgVec <- 3:31
  FVec <- c(0, 1, 1.5)
  
  for (Fmulti in FVec) { # loop over fishing mortality
    Fm <- Fmulti * M 
    message('FM = ', Fm)
    set.seed(101)
    annualF <- Ftrend(1950, 2019, Fm, 'stable', Fcv=0.1, plot=FALSE)
    
    for (binwidth in bwVec) { # loop over bin widths
      message('binwidth = ', binwidth)
      # Simulate population with high-resolution age-structured GTG model 
      SimPop <- GTGpopsim(Linf, K, t0, M, L50, L95, LFS, L5, Vmaxlen, sigmaR,
                          steepness, annualF,alpha, beta, LinfCV, ngtg=1001, maxsd,
                          binwidth)
      # Calculate prob l|a for last simulated year
      DF <- SimPop$df %>% filter(Yr==max(Yr)) # last year
      DF$NVuln <- DF$N * DF$Select
      LenBins <- SimPop$LenBins
      LenMids <- SimPop$LenMids
      Nbins <- length(LenMids)
      ProbLA <- rep(0, Nbins)
      for (l in seq_along(LenMids)) {
        ind <- DF$Length >= LenBins[l] & DF$Length < LenBins[l+1]
        ProbLA[l] <-sum(DF$NVuln[ind])
      }
      ProbLA <- ProbLA/sum(ProbLA) # prob l|a of the vulnerable fished population
      AgeComp <- DF %>% group_by(Age) %>% summarise(vN=sum(N)) # age composition of vulnerable population
      
      for (ngtg in ngtgVec) { # loop over n sub-groups in approx model
        message('ngtg = ', ngtg)
        tt =gc()
        St <- Sys.time()
        SimApprox <- LeesApprox(FVec=annualF, ngtg=ngtg, maxsd, binwidth=binwidth, M,
                                Linf, K, t0,  LFS, L5, Vmaxlen, LinfCV, maxage)
        Elapse <- Sys.time() - St
        probLA_pop <- SimApprox[[1]]
        select_at_length <- SimApprox[[4]]
        
        ProbLA_approx <- rep(0, length(Nbins))
        for (l in seq_along(LenMids)) {
          ProbLA_approx[l] <- sum(probLA_pop[,l] * select_at_length[l]  * AgeComp$vN)
        }
        ProbLA_approx <- ProbLA_approx/sum(ProbLA_approx)
        
        # Calculate Kullback–Leibler divergence between two distributions
        dist <- ProbLA * log(ProbLA/ProbLA_approx)
        dist[!is.finite(dist)] <- 0
        dist <- sum(dist)
        
        xx <- xx + 1
        DFstore[[xx]] <- data.frame(Name=StockNames[st],
                                    Fm=Fm,
                                    Dist=dist,
                                    binwidth=binwidth,
                                    ngtg=ngtg,
                                    Elapse=Elapse,
                                    Fmulti=Fmulti)
      }
    }
  }
}
DF <- do.call("rbind", DFstore)
DF <- DF %>% group_by(Name, Fm, binwidth) %>% mutate(RelDist=Dist/min(Dist))
DF$Fmulti <- as.factor(DF$Fmulti)
DF$binwidth <- as.factor(DF$binwidth)

pt.size <- 1
p1 <- ggplot(DF, aes(x=ngtg, y=Dist,
                     linetype=Fmulti,
                     color=binwidth)) +
  facet_grid(~Name, scales="free_y") +
  expand_limits(y=0, x=1) +
  scale_y_continuous(expand=c(0,0)) +
  geom_line() +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        legend.position = 'none') +
  labs(x="Number of sub-groups (G)",
       y='KL Divergence',
       color="Length class interval",
       # linetype="Length class interval",
       linetype="Relative fishing mortality (F/M)")

maxY <- round(DF$Elapse %>% as.numeric() %>% max(),2)
p2 <- ggplot(DF, aes(x=ngtg, y=Elapse)) +
  facet_grid(~Name, scales="free_y") +
  expand_limits(y=c(0, maxY), x=1) +
  scale_y_continuous(expand=c(0,0)) +
  # geom_point(aes(shape=Fmulti,
  #                color=binwidth)
  #           ) +
  geom_smooth(aes(linetype=Fmulti, color=binwidth),
              method='lm', se=FALSE) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_blank()
  ) +
  labs(x="Number of sub-groups (G)",
       y='Time (seconds)',
       color="Length class interval",
       # linetype="Length class interval",
       linetype="Relative fishing mortality (F/M)")

legend <- cowplot::get_legend(
  # create some space to the left of the legend
  p2 + theme(
    legend.title = element_text(size=8),
    legend.text = element_text(size=8))
)

p2 <- p2 + theme(legend.position = 'none')

Figout <- cowplot::plot_grid(p1, legend, p2, nrow=2, ncol=2, 
                             rel_widths = c(1,0.4),
                             axis='l',
                             align='v',
                             labels=c('a)', '', 'b)'),
                             label_size = 8)
ggsave('Figures/Figure4.png', Figout,
       units='in', height=6, width=8)

