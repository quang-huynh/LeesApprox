
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
Fig4DF <- Figure4()




# ---- Assessment Model (Figure 5) ----

AgeSampSize <- 250
LengthSampSize <- 250

Cobs <- rlnorm(length(annualF), 0, 0.1)
Iobs <- rlnorm(length(annualF), 0, 0.1)


SimPop <- GTGpopsim(Linf, K, t0, M, L50, L95, LFS, L5, Vmaxlen, sigmaR,
                    steepness, annualF,alpha, beta, LinfCV, ngtg=1001, maxsd,
                    binwidth, R0=1E5)

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

AgeSamps <- sapply(2:length(annualF), function(i) 
  rmultinom(n=1, size=AgeSampSize, VAge_Comp[[i+1]]))

AgeSamps <- cbind(rep(0, max(DF$Age)), AgeSamps)
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

# TO DO - Add Obs Error 
Data <- new("Data")
Data@CAL_bins <- SimPop$LenBins
Data@CAL_mids <- SimPop$LenMids
Data@CAL <- array(CAL, dim=c(1, length(annualF), length(SimPop$LenMids)))
Data@CAA <- array(CAA, dim=c(1, length(annualF), max(DF$Age)))
Data@Year <- unique(DF$Yr)
Data@Cat <- matrix(TSData$Catch * Cobs, nrow=1)
Data@Ind <- matrix(TSData$Index * Iobs, nrow=1)

Data@MaxAge <- max(DF$Age)
Data@Mort <- M
Data@vbt0 <- t0
Data@vbK <- K
Data@vbLinf <- Linf
Data@L50 <- L50 
Data@L95 <- L95 
Data@wla <- alpha
Data@CV_vbLinf <- LinfCV
Data@wlb <- beta
Data@steep <- steepness
Data@sigmaR <- sigmaR


CAA_multiplier <- 50
CAL_multiplier <- 0

# Fit Assessment Models 
ngtg_assess <- 3
Mod1 <- SCA_GTG(Data = Data, ngtg=ngtg_assess, 
                CAA_multiplier=CAA_multiplier,
                CAL_multiplier= CAL_multiplier)
Mod2 <- SCA_GTG(Data = Data, ngtg=ngtg_assess, use_LeesEffect = FALSE,
                CAA_multiplier=CAA_multiplier,
                CAL_multiplier= CAL_multiplier)
Mod3 <- SCA(Data=Data,
            CAA_multiplier=CAA_multiplier,
            CAL_multiplier= CAL_multiplier)

plot(annualF, type="l", ylim=c(0, max(annualF*1.5)))
lines(Mod1@FMort, col='blue')
lines(Mod2@FMort, col="green")
lines(Mod3@FMort, col="red")



plot(TSData$Index/max(TSData$Index), type="l", ylim=c(0, 1.5))
lines(Mod1@B_B0, col='blue')
lines(Mod2@B_B0, col="green")
lines(Mod3@B_B0, col="red")


MSEtool::compare_models(Mod1, Mod2, label = c("Lee's Effect", "No Effect"))


source("https://raw.githubusercontent.com/quang-huynh/LeesApprox/master/SCA_GTG_markdown.R")
plot(Mod1)

Mod1@info$data$CAA_n
Obs_C_at_age <- Mod1@Obs_C_at_age
C_at_age <- Mod1@C_at_age
info <- Mod1@info
ind_valid <- rowSums(Obs_C_at_age, na.rm = TRUE) > 0
plot_composition(info$Year[ind_valid], Obs_C_at_age[ind_valid, ], C_at_age[ind_valid, ], plot_type = "annual", ages = NULL, N = info$data$CAA_n[ind_valid])


Mod1@FMort 
Mod2@FMort


# TODO - add some obs error to catch and index 


system.time(
  GTG_3 <- SCA_GTG(Data = SimulatedData, truncate_CAL = FALSE, ngtg=11)
)

# Turn off Lee's Effect. Runtime of 1.1 seconds
system.time(
  GTG_3_noLee <- SCA_GTG(Data = SimulatedData, use_LeesEffect = FALSE, ngtg=11)
)

MSEtool::compare_models(GTG_3, GTG_3_noLee, label = c("Lee's Effect", "No Effect"))
GTG_3@FMort
GTG_3_noLee@FMort

