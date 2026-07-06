## code to prepare `Cond_BaseCase` dataset goes here

rm(list = ls())

dir <- "G:/Shared drives/BM shared/1. Projects/TOF/IOTC-ALB"
load(file.path(dir, 'OM_BaseCase.rda'))

Cond_BaseCase <- list()

# MCMC posterior samples
Cond_BaseCase$mcmcvars         <- mcmcvars

# Biological parameters (one value per MCMC iteration)
Cond_BaseCase$steepness         <- hh
Cond_BaseCase$natural_mortality <- M
Cond_BaseCase$R0                <- R0
Cond_BaseCase$sigmaR            <- sigmar

# Biological schedules (arrays: annual_age × season × sex)
Cond_BaseCase$maturity          <- mata
Cond_BaseCase$weight            <- wta
Cond_BaseCase$ALK               <- pla

# Observed data
Cond_BaseCase$Index             <- I
Cond_BaseCase$LengthComp        <- LF

usethis::use_data(Cond_BaseCase, overwrite = TRUE)
