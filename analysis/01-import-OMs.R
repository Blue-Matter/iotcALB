source('analysis/00-specifications.R')

# Import OMs and Simulate historical fishery dynamics
#
# The six fleets of the ABC OMs are combined into a single fleet with
# MSEtool::CombineFleets() before simulating, to reduce run time and the size
# of the Hist objects.
# The combined fleet reproduces the total fishing mortality-at-age of the six
# fleets. The LL CPUE indices are moved to the Survey slot, each with the
# selectivity of its fleet before combining, so the simulated indices are
# unchanged.

# Saved objects:
#   objects/OM/<OMName>.om                - OM with the six fleets
#   objects/Hist/<OMName>.hist            - Hist with combined fleets (used for the MSE)
#   objects/Hist/<OMName>_all_fleets.hist - Hist with the six fleets (only if AllFleets = TRUE)
#

OM_DF <- data.frame(OM        = OMNames,
                    Import    = TRUE,
                    Simulate  = TRUE,
                    AllFleets = FALSE) # also simulate the six-fleet OM (~1.5 GB Hist)

for (i in seq_len(nrow(OM_DF))) {

  OMName <- OM_DF$OM[i]
  cli::cli_h1('{OMName}')

  OMFile <- file.path('objects/OM', paste0(OMName, '.om'))

  # ---- Import ----
  if (OM_DF$Import[i]) {
    OM <- ImportOM(OMName,
                   Interval      = OMSpecs$Interval,
                   Seasons       = OMSpecs$Seasons,
                   DataLag       = OMSpecs$DataLag,
                   CurrentYear   = OMSpecs$CurrentYear,
                   nYear         = OMSpecs$nYear,
                   pYear         = OMSpecs$pYear,
                   MPStartYear   = OMSpecs$MPStartYear,
                   InterimAdvice = OMSpecs$InterimAdvice,
                   FleetNames    = OMSpecs$FleetNames)

    MSEtool::Save(OM, path = OMFile, overwrite = TRUE)
  }

  # ---- Simulate ----
  if (!OM_DF$Simulate[i]) next

  Versions <- if (OM_DF$AllFleets[i]) c('Combined', 'AllFleets') else 'Combined'

  for (Version in Versions) {
    OM <- readRDS(OMFile)

    # Combine all fleets into one to reduce run time
    if (Version == 'Combined')
      OM <- MSEtool::CombineFleets(OM)

    Hist <- Simulate(OM,
                     control = SimControl(RefPoints       = FALSE,
                                          MGT             = FALSE,
                                          CalcCatchAtSize = FALSE))

    # over-write data used for conditioning with the new updated data
    Hist <- UpdateData(Hist)

    if (Version == 'Combined') {
      ValidateOM(Hist)
      HistFile <- file.path('objects/Hist', paste0(OMName, '.hist'))
    } else {
      ValidateOM(Hist, outdir = 'figures/diagnostics/OM_all_fleets')
      HistFile <- file.path('objects/Hist', paste0(OMName, '_all_fleets.hist'))
    }

    MSEtool::Save(Hist, path = HistFile, overwrite = TRUE)
    rm(Hist); gc()
  }
}
