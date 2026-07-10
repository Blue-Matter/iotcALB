library(MSEtool)
albMSErepo <- '../albMSE'
save_vars  <- ls()

# ---- Fixed Biological Schedules ----

data_path <- file.path(albMSErepo, 'boot/data/alb_abcdata.rda')
data_env  <- new.env()
load(data_path, envir = data_env)

albMSE_Biology <- list(
  ages    = data_env$ages,
  mula    = data_env$mula,
  sdla    = data_env$sdla,
  mulbins = data_env$mulbins,
  wta     = data_env$wta,
  mata    = data_env$mata
)

usethis::use_data(albMSE_Biology, overwrite = TRUE)

Seasons   <- 4L
InitialYr <- 2000L
TermYr    <- 2020L
nYear     <- length(InitialYr:TermYr)
pYear     <- 0L
HistYears <- MSEtool::CalcYears(nYear, pYear, TermYr, Seasons)
nTS       <- nYear * Seasons

StockNames <- c('Female', 'Male')
StockObjectList  <- MSEtool::MakeNamedList(StockNames)

for (st in seq_along(StockObjectList)) {

  stock <- MSEtool::Stock(Name = StockNames[st])

  ## ---- Ages  ----
  annual_age_classes  <- albMSE_Biology$ages
  rep_age_classes     <- rep(annual_age_classes, each = Seasons)
  seas_age_clases     <- seq(min(annual_age_classes),
                             by = 1/Seasons,
                             length.out = length(rep_age_classes))

  MSEtool::Ages(stock) <- MSEtool::Ages(MinAge    = min(seas_age_clases),
                                        MaxAge    = max(seas_age_clases)*Seasons,
                                        Units     = MSEtool::CalcTSUnits(Seasons),
                                        PlusGroup = FALSE)

  AgeClasses <- MSEtool::Classes(stock)
  nAge       <- MSEtool::nAge(stock)
  nAnnualAge <- length(annual_age_classes)

  flat_age_idx  <- c(1L, rep(seq(2L, nAnnualAge), each = 4L), rep(nAnnualAge, 3L))
  flat_seas_idx <- c(4L, rep(1:4, nAnnualAge - 1L), 1:3)

  ## ---- Length -----
  length_at_age <- do.call(cbind, lapply(seq_len(Seasons), \(s)
                                         albMSE_Biology$mula[flat_age_idx, s, st]
  ))

  sd_len_at_age <- do.call(cbind, lapply(seq_len(Seasons), \(s)
                                          albMSE_Biology$sdla[flat_age_idx, s, st]
  ))

  cv_len_at_age <- sd_len_at_age/length_at_age  # seems very small
  midpoints <- albMSE_Biology$mulbins
  by <- diff(midpoints)[1]
  Classes <- seq(midpoints[1] - 0.5*by, by = by, length.out = length(midpoints))

  MSEtool::Length(stock) <- MSEtool::Length(
    MeanAtAge = array(length_at_age,
                      dim = c(1, dim(length_at_age)),
                      dimnames = list(
                        Sim  = 1,
                        Age  = AgeClasses,
                        Year = HistYears[seq_len(ncol(length_at_age))])),
    CVatAge = array(cv_len_at_age,
                      dim = c(1, dim(cv_len_at_age)),
                      dimnames = list(
                        Sim  = 1,
                        Age  = AgeClasses,
                        Year = HistYears[seq_len(ncol(cv_len_at_age))])),
    Classes = Classes,
    Units   = 'cm',
    TruncSD = Inf
    )

  ## ---- Weight ----
  weight <- do.call(cbind, lapply(seq_len(Seasons), \(s)
                                albMSE_Biology$wta[flat_age_idx, s, st]
  ))

  MSEtool::Weight(stock) <- MSEtool::Weight(MeanAtAge = array(weight,
                                                              dim = c(1, dim(weight)),
                                                              dimnames = list(
                                                                Sim  = 1,
                                                                Age  = AgeClasses,
                                                                Year = HistYears[seq_len(ncol(weight))])),
                                            Units = 'kg'
                                            )

  ## ---- Natural Mortality ----
  # from mcmc samples

  ## ---- Maturity ----
  maturity <- do.call(cbind, lapply(seq_len(Seasons), \(s)
                                    albMSE_Biology$mata[flat_age_idx, s, st]
  ))
  MSEtool::Maturity(stock) <- MSEtool::Maturity(MeanAtAge = array(maturity,
                                                                  dim = c(1, dim(maturity)),
                                                                  dimnames = list(
                                                                    Sim  = 1,
                                                                    Age  = AgeClasses,
                                                                    Year = HistYears[seq_len(ncol(maturity))])))

  ## ---- Fecundity ----
  spawning_season <- 3L
  fecundity       <- maturity * weight
  fecundity[,-spawning_season] <- 0  # only produce eggs in `spawning_season`
  MSEtool::Fecundity(stock) <- MSEtool::Fecundity(MeanAtAge = array(fecundity,
                                                                    dim = c(1, dim(fecundity)),
                                                                    dimnames = list(
                                                                      Sim  = 1,
                                                                      Age  = AgeClasses,
                                                                      Year = HistYears[seq_len(ncol(fecundity))])))
  ## ---- Stock-Recruit ----
  # from mcmc samples

  StockObjectList[[st]] <- stock

}

usethis::use_data(StockObjectList, overwrite = TRUE)

rm(list = setdiff(ls(), save_vars))


