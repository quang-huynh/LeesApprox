
# --- Make Figures ----
# Figure1()

# Figure2()

# Figure3()


# --- Simulation Tests ----

SimulateLengths <- function(Linf=100, K=0.1333, t0=0, M=0.2, L50=66, L95=70,
                            L5=40, LFS=50, Vmaxlen=1, sigmaR=0.2, steepness=0.99,
                            annualF, alpha=1E-5, beta=3, LinfCV=0.1,
                            ngtg=1001, maxsd=2, binwidth=5, SampSize=1000) {

  SimPop <- GTGpopsim(Linf, K, t0, M, L50, L95, LFS, L5, Vmaxlen, sigmaR,
                      steepness, annualF,alpha, beta, LinfCV, ngtg, maxsd,
                      binwidth)

  DF <- SimPop$df %>% filter(Yr==max(Yr)) # last year
  DF$NVuln <- DF$N * DF$Select
  LenBins <- SimPop$LenBins
  LenMids <- SimPop$LenMids


  FLengths <- base::sample(DF$Length, size=SampSize, replace=TRUE, prob=DF$NVuln)
  UnFLengths <- base::sample(DF$Length, size=SampSize, replace=TRUE, prob=DF$N)
  out <- list(LenMids=LenMids,
              LenBins=LenBins,
              FishedFreq=cut(FLengths, LenBins) %>% table() %>% as.vector(),
              UnfishedFreq=cut(UnFLengths, LenBins) %>% table() %>% as.vector(),
              DF=DF)
  out
}


library(DLMtool)
DLMextra()
library(DLMextra)
library(dplyr)
library(ggplot2)

loglik <- function(x, p) sum( x * log(p) )

Stocks <- list(Queen_Triggerfish_STT_NOAA,
            Stoplight_Parrotfish_STX_NOAA,
            Yellowtail_Snapper_PR_NOAA)
saveF <- list()

ngtgVec <- c(seq(3, 51, by=2), 1001)
binWidthvec <-  c(2, 5, 10)

for (st in 1:length(Stocks)) {
  if (st == 1) {
    count <- 0
    store <- storeLenComp <- list()
  }
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
  SampSize <- 1000

  annualF <- Ftrend(1950, 2019, 1.5*M, 'stable', Fcv=0.01, plot=FALSE)
  saveF[[st]] <- annualF

  for (bw in binWidthvec) {
    count <- count + 1
    # Simulate length data with high resolution GTG model
    LengthData <- SimulateLengths(Linf, K, t0, M, L50, L95, L5, LFS, Vmaxlen,
                                  sigmaR, steepness, annualF, alpha, beta, LinfCV,
                                  ngtg=1001, maxsd, binwidth=bw, SampSize)

    LenMids <- LengthData$LenMids
    AgeComp <- LengthData$DF %>% group_by(Age) %>% summarise(vN=sum(N))
    Prob <- LengthData$FishedFreq/sum(LengthData$FishedFreq)
    Prob[Prob==0] <- 1E-6

    # NLLhighres <- loglik(LengthData$FishedFreq, Prob)

    nout <- length(ngtgVec)
    Elapse <- NLL <- rep(NA, nout)
    storeLens <- list()
    for (x in seq_along(ngtgVec)) {
      # Simulate length data with Lee's approx method
      St <- Sys.time()
      Sim <- LeesApprox(FVec=annualF, ngtg=ngtgVec[x], maxsd, binwidth=bw, M,
                        Linf, K, t0,  LFS, L5, Vmaxlen, LinfCV, maxage)
      Elapse[x] <- Sys.time() - St

      # Generate length Comp given age comp
      LenComp <- rep(NA, length(LenMids))
      for (l in 1:length(LenMids)) {
        LenComp[l] <- sum(Sim[[1]][,l] * Sim[[4]][l]  * AgeComp$vN)
      }
      LenComp <- LenComp/sum(LenComp) * SampSize

      NLL[x] <- loglik(LenComp, Prob)#/NLLhighres
      storeLens[[x]] <- LenComp

    }
    storeLenComp[[count]] <- data.frame(Length=unlist(storeLens), LenMids=LenMids,
               Ngtg=rep(ngtgVec, each=length(LenMids)),
               binWidth=bw, Stock=st)

    # Elapse <- Elapse/Elapse[length(Elapse)]
    store[[count]] <- data.frame(Elapse=Elapse, binWidth=bw, Stock=st,
                                 Ngtg=ngtgVec, NLL=NLL)
  }
  message(st, '/', length(Stocks))
}
DF <- do.call('rbind', store)

# Relative time of low-res approx compared to high-res model
Max <- DF %>% group_by(Stock, binWidth) %>% filter(Ngtg==1001) %>%
  summarise(Maxtime=Elapse, MaxNLL=NLL)
DF <- left_join(DF, Max, by=c("Stock", "binWidth"))

Elapse <- DF %>% group_by(Stock, binWidth, Ngtg) %>%
  summarise(RelElapse=Elapse/Maxtime)

Elapse %>% group_by(Stock, binWidth) %>% summarise(mean(RelElapse))


Elapse %>% filter(Ngtg  != 1001) %>% group_by(Stock) %>%
  summarize(mean=1-mean(RelElapse))


dat <- Elapse %>% group_by(Stock, Ngtg) %>%
  filter(Ngtg  != 1001) %>% summarize(mean=mean(RelElapse))

dat %>% filter(Stock==1) %>% mutate(1/mean) %>% data.frame()
dat %>% filter(Stock==2) %>% mutate(1/mean) %>% data.frame()
dat %>% filter(Stock==3) %>% mutate(1/mean) %>% data.frame()


# Compare NLL
NLLdat <- DF %>% group_by(Stock, Ngtg, binWidth) %>%  filter(Ngtg!=1001) %>%
  summarise(relNLL=-NLL/-MaxNLL) %>% ungroup() %>%
  group_by(Stock, Ngtg) %>% summarise(mean=mean(relNLL))

Species <- data.frame(Stock=1:3, Species=c("Yellowtail snapper", "Queen triggerfish", "Stoplight parrotfish"))
NLLdat <- left_join(NLLdat, Species, by='Stock')

P1 <- ggplot(NLLdat, aes(x=Ngtg, y=mean, linetype=Species, color=Species)) + geom_line() +
  geom_point(size=2) +
  theme_classic() +
  labs(x="Number of sub-groups", y="Relative Log-likelihood")

ggsave("Figures/Figure4.png", P1, width=5, height=3)





# Plot Size Comps
cols <- c("darkblue", "darkgray")
DF2 <- do.call('rbind', storeLenComp) %>% left_join(., Species)

plotDat <- DF2 %>% filter(binWidth==10, Ngtg%in% c(3, 5, 11, 51, 1001))

pdat <- plotDat %>% filter(Stock==1)
pdat2 <- pdat %>% filter(Ngtg ==1001)
pdat3 <- pdat %>% filter(Ngtg !=1001)
pdat3$NGTG <- pdat3$Ngtg
pdat3$Model <- "Approximation"
pdat2a <- rbind(pdat2, pdat2, pdat2, pdat2)
pdat2a$Model <- "High Resolution"
pdat2a$NGTG <- rep(unique(pdat3$Ngtg), each=nrow(pdat2))
pdat4 <- rbind(pdat3, pdat2a)
pdat4 <- pdat4 %>% filter(Length >=1)
P2a <- ggplot(pdat4, aes(x=LenMids, y=Length, fill=Model)) + facet_grid(Species~NGTG) +
  geom_bar(stat="identity",position = "identity", alpha=.4) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  labs(x="", y="Frequency")

pdat <- plotDat %>% filter(Stock==2)
pdat2 <- pdat %>% filter(Ngtg ==1001)
pdat3 <- pdat %>% filter(Ngtg !=1001)
pdat3$NGTG <- pdat3$Ngtg
pdat3$Model <- "Approximation"
pdat2a <- rbind(pdat2, pdat2, pdat2, pdat2)
pdat2a$Model <- "High Resolution"
pdat2a$NGTG <- rep(unique(pdat3$Ngtg), each=nrow(pdat2))
pdat4 <- rbind(pdat3, pdat2a)
pdat4 <- pdat4 %>% filter(Length >=1)
P2b <- ggplot(pdat4, aes(x=LenMids, y=Length, fill=Model)) + facet_grid(Species~NGTG) +
  geom_bar(stat="identity",position = "identity", alpha=.4) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(x="", y="Frequency")

pdat <- plotDat %>% filter(Stock==3)
pdat2 <- pdat %>% filter(Ngtg ==1001)
pdat3 <- pdat %>% filter(Ngtg !=1001)
pdat3$NGTG <- pdat3$Ngtg
pdat3$Model <- "Approximation"
pdat2a <- rbind(pdat2, pdat2, pdat2, pdat2)
pdat2a$Model <- "High Resolution"
pdat2a$NGTG <- rep(unique(pdat3$Ngtg), each=nrow(pdat2))
pdat4 <- rbind(pdat3, pdat2a)
pdat4 <- pdat4 %>% filter(Length >=1)
P2c <- ggplot(pdat4, aes(x=LenMids, y=Length, fill=Model)) + facet_grid(Species~NGTG) +
  geom_bar(stat="identity",position = "identity", alpha=.4) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(x="Length", y="Frequency")

P2 <- DLMtool::join_plots(list(P2a, P2b, P2c), nrow=3, ncol=1, position = "bottom")


ggsave("Figures/Figure5.png", P2, width=9, height=6)






