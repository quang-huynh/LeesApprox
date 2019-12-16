
# --- Code for Lee's Approximation Manuscript ----
# A. Hordyk
# December 2019

library(LeesApproxTMB)
library(dplyr)
library(Rcpp)
library(ggplot2)
library(cowplot)

DLMextra()
library(DLMextra)

setwd('Manuscript') 
source('Functions.r')
Rcpp::sourceCpp('src/LeesApprox.cpp') # load CPP version of Lee's approx

# ---- Figure 1 ----

source('Figure_1.r')

Figure1()

# ---- Figure 2 ----
pars <- list()
pars$ngtg <- 7
pars$LinfCV <- 0.1
pars$Linf <- 100
pars$maxsd <- 2
pars$binwidth <- 5

pars$K <- 0.25
pars$t0 <- 0
pars$M <- 0.2
pars$maxage <- ceiling(-log(0.01)/pars$M)
pars$Ages <- 1:pars$maxage

pars$L5 <- 30
pars$LFS <- 60
pars$Vmax <- 0.6

pars$yr.st <- 1950
pars$yr.end <- 2017
pars$FVec <- c(0, 0.2, 0.6)


source('Figure_2.r')

Figure2(pars)

# ---- Figure 3 ----
source('Figure_3.r')

Figure3(pars)


# ---- Sensitivity Tests (Figure 4) ----
source('Figure_4.r')

Figure4()


# ---- Assessment Model (Figure 5) ----

SimPop <- GTGpopsim(Linf, K, t0, M, L50, L95, LFS, L5, Vmaxlen, sigmaR,
                    steepness, annualF,alpha, beta, LinfCV, ngtg=1001, maxsd,
                    binwidth, R0=1E5)

DF <- SimPop[[1]]

DF <- DF %>% group_by(Yr) %>% mutate(Catch=sum(CAA*Weight),
                                     TotalB=sum(N*Weight))


Catch <- DF %>% select(Yr, Catch) %>% distinct()
Index <- DF %>% select(Yr, TotalB) %>% distinct()

library(purrr)
library(tidyr)
size <- 100
DF$vulnN <- DF$N * DF$Select
AnnualSamples <- DF %>% group_by(Yr) %>% 
  tidyr::nest() %>%
  ungroup() %>%
  mutate(nsamp=size) %>%
  mutate(rows=purrr::map2(data, nsamp, sample_n,
                         replace=TRUE, 
                         weight=vulnN)) %>%
  select(-data) %>%
  unnest(rows)

dim(AnnualSamples)
head(AnnualSamples)

CAA_samples <- AnnualSamples %>% group_by(Yr, Age) %>% summarize(CAA=sum(CAA))

tt <- CAA_samples %>% tidyr::pivot_wider(names_from='Yr',
                                   values_from="CAA")

CAA_samples$Yr %>% unique()
CAA_samples$Age %>% unique()


CAA_samples %>% filter(Yr==1)


head(tt)
dim(tt)


# Catch at age in numbers?

plot(Catch, type='b')
plot(Index, type='b')

head(DF)

Data <- new("Data")
Data@CAL_bins <- SimPop$LenBins
Data@CAL_mids <- SimPop$LenMids
Data@Year <- unique(DF$Yr)
Data@Cat <- matrix(Catch$Catch, nrow=1)
Data@Ind <- matrix(Index$Index, nrow=1)

Data@vbt0 <- t0
Data@vbK <- K
Data@vbLinf <- Linf
Data@L50 <- L50 
Data@L95 <- L95 


# TODO - add some obs error to catch and index 

# Catch data

# Index 

# CAA 

# CAL 

# Life-history






?SCA_GTG()



 
# Generate fishery data to test assessment







# for (bw in binWidthvec) {
  # count <- count + 1
  # Simulate length data with high resolution GTG model
LengthData <- SimulateLengths(Linf, K, t0, M, L50, L95, L5, LFS, Vmaxlen,
                              sigmaR, steepness, annualF, alpha, beta, LinfCV,
                              ngtg=1001, maxsd, binwidth=bw, SampSize)
  




  probL <- rep(0, length(LenMids))
  for (l in seq_along(LenMids)) {
    ind <- DF$Length >= LenBins[l] & DF$Length < LenBins[l+1]
    probL[l] <- sum(DF$NVuln[ind])
  }

  probL <- probL/sum(probL)



  out <- list(LenMids=LenMids,
              LenBins=LenBins,
              FishedFreq=probL,
              DF=DF)
  out
}




# loglik <- function(x, p) sum( x * log(p) )


saveF <- list()
Stock_Names <- c("Queen triggerfish", "Stoplight parrotfish", "Yellowtail snapper")

ngtgVec <- c(seq(3, 51, by=2), 1001)
binWidthvec <-  c(2, 5, 10)

for (st in 1:length(Stocks)) {
  if (st == 1) {
    count <- 0
    store <- storeLenComp <- list()
  }
  Stock <- Stocks[[st]]
  Name <- Stock_Names[st]
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
  SampSize <- 1000
  maxsd <- 2


  annualF <- Ftrend(1950, 2019, 1.5*M, 'stable', Fcv=0.01, plot=FALSE)
  saveF[[st]] <- annualF

  for (bw in binWidthvec) {
    count <- count + 1
    # Simulate length data with high resolution GTG model
    LengthData <- SimulateLengths(Linf, K, t0, M, L50, L95, L5, LFS, Vmaxlen,
                                  sigmaR, steepness, annualF, alpha, beta, LinfCV,
                                  ngtg=1001, maxsd, binwidth=bw)

    LenMids <- LengthData$LenMids
    AgeComp <- LengthData$DF %>% group_by(Age) %>% summarise(vN=sum(N))
    # AgeComp <- LengthData$FAges
    Prob <- LengthData$FishedFreq/sum(LengthData$FishedFreq)
    Prob[Prob==0] <- 1E-6
    Prob <- Prob/sum(Prob)

    nout <- length(ngtgVec)
    NLL <- rep(NA, nout)
    Elapse <- rep(NA, nout)
    storeLens <- list()
    for (x in seq_along(ngtgVec)) {
      # Simulate length data with Lee's approx method
      gc()
      St <- Sys.time()
      Sim <- LeesApprox(FVec=annualF, ngtg=ngtgVec[x], maxsd, binwidth=bw, M,
                        Linf, K, t0,  LFS, L5, Vmaxlen, LinfCV, maxage)
      Elapse[x] <- Sys.time() - St

      # Generate length Comp given age comp
      LenComp <- rep(NA, length(LenMids))
      for (l in 1:length(LenMids)) {
        LenComp[l] <- sum(Sim[[1]][,l] * Sim[[4]][l]  * AgeComp$vN)
      }
      # LenComp <- LenComp/sum(LenComp) * SampSize

      LenComp <- LenComp/sum(LenComp)
      LenComp[LenComp==0] <- 1E-6

      NLL[x] <- sum(Prob * log(LenComp/Prob))

      storeLens[[x]] <- LenComp * SampSize

    }
    storeLenComp[[count]] <- data.frame(Length=unlist(storeLens), LenMids=LenMids,
               Ngtg=rep(ngtgVec, each=length(LenMids)),
               binWidth=bw, Stock=st)

    # Elapse <- Elapse/Elapse[length(Elapse)]
    store[[count]] <- data.frame(Name=Name, Elapse=Elapse, binWidth=bw, Stock=st,
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
  mutate(RelElapse=Elapse/Maxtime)

Elapse <- Elapse %>% filter(Ngtg  != 1001)
ggplot(Elapse, aes(x=Ngtg, y=RelElapse, linetype=as.factor(binWidth))) +
  facet_wrap(~Stock) +
  geom_line() + geom_point() +
  geom_smooth(method='lm', formula= y~x, se=FALSE) +
  expand_limits(y=0) + theme_classic()


mod <- lm(Elapse ~ Ngtg + binWidth, data=Elapse)

n <- 3
bw <- 10
coef(mod)[1] + n*coef(mod)[2] + bw*coef(mod)[3]

predict(mod, newdata=data.frame(Ngtg=n, binWidth=bw))



Elapse %>% group_by(Stock, binWidth) %>% summarise(mean(RelElapse))


Elapse %>% filter(Ngtg  != 1001) %>% group_by(Stock) %>%
  summarize(mean=1-mean(RelElapse))


dat <- Elapse %>% group_by(Stock, Ngtg) %>%
  filter(Ngtg  != 1001) %>% summarize(mean=mean(RelElapse))

dat %>% filter(Stock==1) %>% mutate(1/mean) %>% data.frame()
dat %>% filter(Stock==2) %>% mutate(1/mean) %>% data.frame()
dat %>% filter(Stock==3) %>% mutate(1/mean) %>% data.frame()


# Compare NLL
NLLdat <- DF %>% group_by(Name, Ngtg, binWidth) %>%  filter(Ngtg!=1001) %>%
  summarise(relNLL=-NLL/-MaxNLL) %>% ungroup()

ggplot(NLLdat, aes(x=Ngtg, y=relNLL, color=as.factor(binWidth))) +
  facet_grid(~Name) +
  geom_line()

DF %>% filter(Name=='Queen triggerfish', binWidth==10)






NLLdat <- DF %>% group_by(Name, Ngtg, binWidth) %>%  filter(Ngtg!=1001) %>%
  summarise(relNLL=-NLL/-MaxNLL) %>% ungroup() %>%
  group_by(Name, Ngtg) %>% summarise(mean=mean(relNLL))

P1 <- ggplot(NLLdat, aes(x=Ngtg, y=mean, linetype=Name, color=Name)) + geom_line() +
  geom_point(size=2) +
  theme_classic() +
  labs(x="Number of sub-groups", y="Relative Log-likelihood")

ggsave("Figures/Figure4.png", P1, width=5, height=3)





# Plot Size Comps
cols <- c("black", "darkgray")
DF2 <- do.call('rbind', storeLenComp) %>% left_join(., Species)

alpha <- 0.6
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
  geom_bar(stat="identity",position = "identity", alpha=alpha) +
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
  geom_bar(stat="identity",position = "identity", alpha=alpha) +
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
  geom_bar(stat="identity",position = "identity", alpha=alpha) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(x="Length", y="Frequency")

P2 <- DLMtool::join_plots(list(P2a, P2b, P2c), nrow=3, ncol=1, position = "bottom")


ggsave("Figures/Figure5.png", P2, width=9, height=6)


# --- Fit to historical data ----

Linf=100; K=0.1333; t0=0; M=0.2; L50=66; L95=70;
L5=40; LFS=50; Vmaxlen=1; sigmaR=0.2; steepness=0.99;
alpha=1E-5; beta=3; LinfCV=0.1;
ngtg=1001; maxsd=2; binwidth=5; SampSize=1000
#
library(dplyr)
LengthData <- SimulateLengths(Linf, K, t0, M, L50, L95, L5, LFS, Vmaxlen,
                              sigmaR, steepness, annualF, alpha, beta, LinfCV,
                              ngtg=1001, maxsd, binwidth=binwidth)


Data <- new("Data")

Data@

GTG_1 <- SCA_GTG(Data = Data, ngtg=15)

