library(iotcALB)

OMNames <- c('OM5a', 'OM5b', 'OM6b')

for (OMName in OMNames) {

  cli::cli_h1('{OMName}')

  OM <- ImportOM(OMName)
  MSEtool::Save(OM, path = file.path('objects/OM', paste0(OMName, '.om')), overwrite = TRUE)

  Hist <- Simulate(OM)
  ValidateOM(Hist)
  MSEtool::Save(Hist, path = file.path('objects/Hist', paste0(OMName, '.hist')), overwrite = TRUE)
}

