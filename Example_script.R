
devtools::install_github("quang-huynh/LeesApprox/LeesApproxTMB")
library(LeesApproxTMB)
library(MSEtool)
?SCA_GTG

# Default with 3 GTG's and use_Lees' effect = TRUE
# Interpolation is always used.
# SCA_GTG fits to both age and length comps.
# By default, no age comps are used, argument: CAA_multiplier = 0.

# This took about 20 seconds to run
system.time(
GTG_3 <- SCA_GTG(Data = SimulatedData)
)

# Turn off Lee's Effect. Runtime of 1.1 seconds
system.time(
GTG_3_noLee <- SCA_GTG(Data = SimulatedData, use_LeesEffect = FALSE)
)

# Compare the model with and without Lee's Effect
# As expected, we estimate lower F's and higher stock sizes with Lee's Effect
MSEtool::compare_models(GTG_3, GTG_3_noLee, label = c("Lee's Effect", "No Effect"))

# Markdown reporting for GTG_3
source("https://raw.githubusercontent.com/quang-huynh/LeesApprox/master/SCA_GTG_markdown.R")
plot(GTG_3)

# Let's increase to 31 GTGs. 2.5 minute runtime
system.time(
GTG_31 <- SCA_GTG(Data = SimulatedData, ngtg = 31)
)

system.time(
  GTG_501 <- SCA_GTG(Data = SimulatedData, ngtg = 501, max_sd_gtg = 3)
)

system.time(
  GTG_501_noLee <- SCA_GTG(Data = SimulatedData, ngtg = 501, use_LeesEffect = FALSE)
)

# Compare the models with 3 and 31 GTGs. Three works pretty good!
compare_models(GTG_3, GTG_31, label = c("3 GTG", "31 GTG"))

# Let's compare with an SCA that uses age comps
system.time(
age_SCA <- SCA(Data = SimulatedData)
)
compare_models(GTG_3, GTG_31, GTG_501, label = c("3 GTG", "31 GTG", "501 GTG"))
