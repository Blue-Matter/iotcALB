library(iotcALB)

# ---- Base Case ----

Name <- 'Base Case'

# Import and Save OM
BaseCase <- ImportOM(object = Cond_BaseCase, Name = Name)
Save(BaseCase, path = paste0('objects/OM/', Name, '.om'), overwrite = TRUE)

# Simulate Historical Fishery
Hist_BaseCase <- Simulate(BaseCase)

# Validate
ValidateOM(Hist_BaseCase, object = Cond_BaseCase)

# Save Hist
Save(Hist_BaseCase, path = paste0('objects/Hist/', Name, '.hist'), overwrite = TRUE)


# ---- Robustness:  ----

# TBD
