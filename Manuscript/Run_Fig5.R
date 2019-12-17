# devtools::install_github('adrianhordyk/LeesApprox')

library(LeesApproxTMB)
library(dplyr)
library(DLMextra)
source('Manuscript/functions.r')
source('Manuscript/Figure_5.r')

CVs <- 0.1 # obs error on index and catches
LH_CV <- 0.05 # obs error on life-history parameters

LengthSampSize <- 100
AgeSampSize <- 100
Years <- -10 # index, CAA and CAL data for last 10 years


Stock <- 2 

genData <- GenerateData(Stock=Stock, Years=Years, 
                        CobCV=CVs, IobCV=CVs,
                        LH_CV=LH_CV,
                        LengthSampSize=LengthSampSize,
                        AgeSampSize=AgeSampSize,
                        Fmulti=1, nbins=30)

annualF <- genData$annualF
Data <- genData$Data
LHpars <- genData$LHpars


vulnerability <- 'logistic' # 'dome'

control <- list(iter.max = 5e+05, eval.max = 7e+05)
# Fit Assessment Models - Only CAA data 
CAA_multiplier <- 50
CAL_multiplier <- 0
ngtg_assess <- 11

# with GTG approx
Mod1 <- SCA_GTG(Data = Data, ngtg=ngtg_assess, 
                CAA_multiplier=CAA_multiplier,
                CAL_multiplier= CAL_multiplier, 
                control=control)
# without GTG approx
Mod2 <- SCA_GTG(Data = Data, ngtg=ngtg_assess, 
                use_LeesEffect = FALSE,
                CAA_multiplier=CAA_multiplier,
                CAL_multiplier= CAL_multiplier,
                control=control)
# MSEtool::SCA
Mod3 <- SCA(Data=Data,
            CAA_multiplier=CAA_multiplier,
            CAL_multiplier= CAL_multiplier,
            control=control)

data.frame(Mod=c('Mod1', 'Mod2', 'Mod3'),
           conv=c(Mod1@conv,Mod2@conv, Mod3@conv))
           
cbind(Mod1@F_FMSY, Mod2@F_FMSY, Mod3@F_FMSY)
cbind(Mod1@opt, Mod2@opt,Mod3@opt)   

plot(annualF, type="l", ylim=c(0, max(annualF*1.5)))
lines(Mod1@FMort, col='blue')
lines(Mod2@FMort, col="green")
lines(Mod3@FMort, col="red")

plot(genData$Index/genData$Index[1], type="l", ylim=c(0, 1.5))
lines(Mod1@B_B0, col='blue')
lines(Mod2@B_B0, col="green")
lines(Mod3@B_B0, col="red")


source("Manuscript/SCA_GTG_markdown.r")
plot(Mod1)
plot(Mod2)




MSEtool::compare_models(Mod2, Mod3)



matplot(t(Mod1@Selectivity), type="l")
matplot(t(Mod2@Selectivity), type="l")
matplot(t(Mod3@Selectivity), type="l")





# Fit Assessment Models - Only CAL data 
CAA_multiplier <- 0
CAL_multiplier <- 50
ngtg_assess <- 11
Mod1 <- SCA_GTG(Data = Data, ngtg=ngtg_assess, 
                CAA_multiplier=CAA_multiplier,
                CAL_multiplier= CAL_multiplier)

Mod2 <- SCA_GTG(Data = Data, ngtg=ngtg_assess, 
                use_LeesEffect = FALSE,
                CAA_multiplier=CAA_multiplier,
                CAL_multiplier= CAL_multiplier)


plot(annualF, type="l", ylim=c(0, max(annualF*1.5)))
lines(Mod1@FMort, col='blue')
lines(Mod2@FMort, col="green")


matplot(t(Mod1@Selectivity), type="l")
matplot(t(Mod2@Selectivity), type="l")















# ---- Loop over Stocks and run nsim estimates ----

xx <- 0 
nsim <- 20
outpars <- list()
for (st in 1:3) {
  for (i in 1:nsim) {
    message(i, '/', nsim)
    genData <- GenerateData(st=st, Years=Years, 
                            CobCV=CVs, IobCV=CVs,
                            LengthSampSize=LengthSampSize,
                            AgeSampSize=AgeSampSize,
                            Fmulti=1, nbins=30)
    
    annualF <- genData$annualF
    Data <- genData$Data
    LHpars <- genData$LHpars
    
    
    vulnerability <- 'logistic' # 'dome'
    # Fit Assessment Models - Only CAL data 
    CAA_multiplier <- 0
    CAL_multiplier <- 50
    ngtg_assess <- 11
    Data@Ind[Data@Ind>0] <- NA
    gc()
    start <- Sys.time()
    Mod1 <- SCA_GTG(Data = Data, ngtg=ngtg_assess, 
                    CAA_multiplier=CAA_multiplier,
                    CAL_multiplier=CAL_multiplier,
                    vulnerability=vulnerability )
    GTG_Elapse <- Sys.time() - start
    
    gc()
    start <- Sys.time()
    Mod2 <- SCA_GTG(Data = Data, ngtg=ngtg_assess, 
                    use_LeesEffect = FALSE,
                    CAA_multiplier=CAA_multiplier,
                    CAL_multiplier= CAL_multiplier,
                    vulnerability=vulnerability )
    nGTG_Elapse <- Sys.time() - start
    
    termF <- annualF[length(annualF)]
    estF_GTG <- Mod1@FMort[length(Mod1@FMort)]
    estF_nGTG <- Mod2@FMort[length(Mod2@FMort)]
    
    termD <- genData$Index[length(genData$Index)]/genData$Index[1]
    estD_GTG <- Mod1@B_B0[length(Mod1@B_B0)]
    estD_nGTG <- Mod2@B_B0[length(Mod2@B_B0)]
    
    xx <- xx + 1 
    outpars[[xx]] <- data.frame(LHpars, termF, estF_GTG, estF_nGTG, 
                                termD, estD_GTG, estD_nGTG, 
                                GTG_Elapse, nGTG_Elapse)
  }
}

Results <- do.call('rbind', outpars)

saveRDS(Results, 'Manuscript/AssessResults.rda')



Results <- Results %>% mutate(restF_GTG=estF_GTG/termF,
                              rest_nGTG=estF_nGTG/termF)
tt <- tidyr::pivot_longer(Results, cols=c('restF_GTG', 'rest_nGTG'))

library(ggplot2)
ggplot(tt, aes(x=name, y=value)) + facet_grid(~Stock) +
  geom_boxplot()


  




plot(annualF, type="l", ylim=c(0, max(annualF*1.5)))
lines(Mod1@FMort, col='blue')
lines(Mod2@FMort, col="green")

cbind(annualF, Mod1@FMort, Mod2@FMort)


matplot(t(Mod1@Selectivity), type="l", ylim=c(0,1))
matplot(t(Mod2@Selectivity), type="l", ylim=c(0,1))


plot(genData$Index/genData$Index[1], type="l", ylim=c(0, 1.5))
lines(Mod1@B_B0, col='blue')
lines(Mod2@B_B0, col="green")


# source("Manuscript/SCA_GTG_markdown.r")
# plot(Mod1)
# plot(Mod2)


# Fit Assessment Models - Only CAA data 
CAA_multiplier <- 50
CAL_multiplier <- 0
ngtg_assess <- 11
Mod1 <- SCA_GTG(Data = Data, ngtg=ngtg_assess, 
                CAA_multiplier=CAA_multiplier,
                CAL_multiplier= CAL_multiplier)

Mod2 <- SCA_GTG(Data = Data, ngtg=ngtg_assess, 
                use_LeesEffect = FALSE,
                CAA_multiplier=CAA_multiplier,
                CAL_multiplier= CAL_multiplier)

Mod3 <- SCA(Data=Data,
            CAA_multiplier=CAA_multiplier,
            CAL_multiplier= CAL_multiplier)

plot(annualF, type="l", ylim=c(0, max(annualF*1.5)))
lines(Mod1@FMort, col='blue')
lines(Mod2@FMort, col="green")
lines(Mod3@FMort, col="red")

matplot(t(Mod1@Selectivity), type="l")
matplot(t(Mod2@Selectivity), type="l")
matplot(t(Mod3@Selectivity), type="l")






