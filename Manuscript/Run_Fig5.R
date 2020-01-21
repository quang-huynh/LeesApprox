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
DatYears <- 10 # NA = data from all years, otherwise index, CAA and CAL data for last `DatYears` years


Stock <- 2

genData <- GenerateData(Stock=Stock, DatYears=DatYears, 
                        CobCV=CVs, IobCV=CVs,
                        LH_CV=LH_CV,
                        LengthSampSize=LengthSampSize,
                        AgeSampSize=AgeSampSize,
                        Fmulti=1, nbins=30)

annualF <- genData$annualF
Data <- genData$Data
LHpars <- genData$LHpars


vulnerability <- 'logistic' # 'dome'
fix_sigma <- TRUE 

control <- list(iter.max = 5e+05, eval.max = 7e+05)
# Fit Assessment Models - Only CAA data 
CAA_multiplier <- 40
CAL_multiplier <- 0
ngtg_assess <- 11

# with GTG approx
Mod1 <- SCA_GTG(Data = Data, ngtg=ngtg_assess, vulnerability = vulnerability,
                CAA_multiplier=CAA_multiplier, rescale = 1,
                CAL_multiplier= CAL_multiplier, start = list(omega = 0.01),
                control=control, 
                fix_sigma =fix_sigma)
# without GTG approx
Mod2 <- SCA_GTG(Data = Data, ngtg=ngtg_assess, vulnerability = vulnerability,
                use_LeesEffect = FALSE,
                CAA_multiplier=CAA_multiplier, rescale = 1,
                CAL_multiplier= CAL_multiplier, start = list(omega = 0.01),
                control=control,
                fix_sigma =fix_sigma)
# MSEtool::SCA
Mod3 <- SCA(Data=Data,
            CAA_multiplier=CAA_multiplier, vulnerability = vulnerability,
            CAL_multiplier= CAL_multiplier, start = list(omega = 0.01), rescale = 1,
            control=control,
            fix_sigma=fix_sigma)


data.frame(Mod=c('Mod1', 'Mod2', 'Mod3'),
           conv=c(Mod1@conv,Mod2@conv, Mod3@conv))
         
# Maximum gradient
sapply(c(Mod1, Mod2, Mod3), function(x) max(abs(x@SD$gradient.fixed)))

# Any zero gradient at starting values
Mod1@obj$gr(Mod1@obj$par)
Mod2@obj$gr(Mod2@obj$par)
Mod3@obj$gr(Mod3@obj$par)

# Any zero gradient at final parameter values
Mod1@obj$gr(Mod1@opt$par)
Mod2@obj$gr(Mod2@opt$par)
Mod3@obj$gr(Mod3@opt$par)

# FMSY
Mod1@FMSY  
Mod2@FMSY
Mod3@FMSY

# Negative log-likelihood components
Mod1@NLL[-4] # The 4th element is for length comps which should be zero. Remove to align with Mod3 vector
Mod2@NLL[-4]
Mod3@NLL

# Correlation between R0, and the 2 selectivity parameters
(Mod1@SD$cov.fixed %>% cov2cor() %>% round(3))[1:3, 1:3]
(Mod2@SD$cov.fixed %>% cov2cor() %>% round(3))[1:3, 1:3]
(Mod3@SD$cov.fixed %>% cov2cor() %>% round(3))[1:3, 1:3]


# F/FMSY
matplot(cbind(Mod1@F_FMSY, Mod2@F_FMSY, Mod3@F_FMSY), type="l")
cbind(Mod1@opt, Mod2@opt,Mod3@opt)   


plot(annualF, type="l", ylim=c(0, max(annualF*1.5)))
lines(Mod1@FMort, col='blue')
lines(Mod2@FMort, col="green") # expect Mod2 and Mod3 to be similiar
lines(Mod3@FMort, col="red")

plot(genData$Index/genData$Index[1], type="l", ylim=c(0, 1.5))
lines(Mod1@B_B0, col='blue')
lines(Mod2@B_B0, col="green") # B_B0 > 1
lines(Mod3@B_B0, col="red") # B_B0 > 1



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


  