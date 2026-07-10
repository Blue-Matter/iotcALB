library(MSEtool)
albMSErepo <- '../albMSE'
model_dir  <- file.path(albMSErepo, 'data')
save_vars  <- ls()

om_dirs <- list.dirs(model_dir, recursive = FALSE)

for (om_dir in om_dirs) {

  om_nm  <- basename(om_dir)
  suffix <- sub('^om', '', om_nm)
  OM_nm  <- paste0('CondData_OM', suffix)

  path <- file.path(om_dir, paste0('mcvars_abc', suffix, '.rda'))

  e <- new.env()
  load(path, envir = e)

  assign(OM_nm, e$mcvars)

  save(list = OM_nm,
       file = file.path('data', paste0(OM_nm, '.rda')),
       compress = 'bzip2')
}

rm(list = setdiff(ls(), save_vars))
