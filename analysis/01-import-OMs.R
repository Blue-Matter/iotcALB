library(iotcALB)

# ---- Base Case ----

Name <- 'Base Case'

BaseCase <- ImportOM(object = Cond_BaseCase, Name = Name)
Save(BaseCase, path = paste0('objects/OM/', Name, '.om'), overwrite = TRUE)

Hist_BaseCase <- Simulate(BaseCase)
ValidateOM(Hist_BaseCase, object = Cond_BaseCase)
Save(Hist_BaseCase, path = paste0('objects/Hist/', Name, '.hist'), overwrite = TRUE)


# ---- Robustness: [Name] ----

