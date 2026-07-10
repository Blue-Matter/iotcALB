library(MSEtool)
albMSErepo <- '../albMSE'
save_vars  <- ls()

Seasons   <- 4L
InitialYr <- 2000L
TermYr    <- 2020L
nYear     <- length(InitialYr:TermYr)
pYear     <- 0L
HistYears <- MSEtool::CalcYears(nYear, pYear, TermYr, Seasons)
nTS       <- nYear * Seasons


nFleet     <- 6
FleetNames <- c(paste0('LL', 1:4), 'PS', 'Other')

data_path  <- file.path(albMSErepo, 'boot/data/alb_abcdata.rda')
data_env   <- new.env()
load(data_path, envir = data_env)

albMSE_Data <- MSEtool::Data(Name = 'IOTC Albacore Data',
                      CommonName = 'Albacore',
                      Species = 'Thunnus alalunga',
                      Years = HistYears,
                      YearLH = max(as.integer(HistYears)),
                      Seasons = Seasons
)


# ---- Catches ----
catch_array <- array(NA,
                     dim = c(nTS, nFleet),
                     dimnames = list(
                       Year = HistYears,
                       Fleet = FleetNames
                     )
)
C_perm <- aperm(data_env$C, c(2, 1, 3))          # [Season, Year, Fleet]
catch_array[] <- matrix(C_perm, nrow = nTS, ncol = nFleet)

MSEtool::Landings(albMSE_Data) <- MSEtool::CatchData(Name  =  FleetNames,
                                              Value = catch_array,
                                              Units = rep('t', nFleet)
                                              )

# ---- Indices ----
ind_data  <- data_env$I
n_ind     <- dim(ind_data)[3]
ind_array <- array(NA,
                     dim = c(nTS, n_ind),
                     dimnames = list(
                       Year = HistYears,
                       Fleet = FleetNames[seq_len(n_ind)]
                     )
)

I_perm <- aperm(ind_data, c(2, 1, 3))          # [Season, Year, Fleet]
ind_array[] <- matrix(I_perm, nrow = nTS, ncol = n_ind)

# standdrdize to mean 1
ind_array <- ind_array/matrix(apply(ind_array, 'Fleet', mean, na.rm=TRUE), nTS, n_ind, byrow=T)

MSEtool::CPUE(albMSE_Data) <- MSEtool::IndicesData(Name  = FleetNames[seq_len(n_ind)],
                                            Value = ind_array)

# ---- Length Frequencies ----
# leaving for now as not used by CMPs


# ---- Save object ----
usethis::use_data(albMSE_Data, overwrite = TRUE)

rm(list = setdiff(ls(), save_vars))
