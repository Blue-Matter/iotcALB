source('analysis/00-specifications.R')

# Import OMs and Simulate historical fishery dynamics


for (i in seq_along(OMNames)) {

  OMName <- OMNames[i]
  cli::cli_h1('{OMName}')

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

  MSEtool::Save(OM, path = file.path('objects/OM', paste0(OMName, '.om')),
                overwrite = TRUE)

  Hist <- Simulate(OM,
                   control = SimControl(RefPoints       = FALSE,
                                        MGT             = FALSE,
                                        CalcCatchAtSize = FALSE))


  Hist <- UpdateData(Hist)
  ValidateOM(Hist)

  MSEtool::Save(Hist, path = file.path('objects/Hist', paste0(OMName, '.hist')), overwrite = TRUE)
}


