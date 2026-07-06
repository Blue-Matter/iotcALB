library(iotcALB)

# ---- Base Case ----

Name <- 'Base Case'

BaseCase <- ImportOM(object = Cond_BaseCase, Name = Name)
Save(BaseCase, path = paste0('objects/OM/', Name, '.om'), overwrite = TRUE)

Hist_BaseCase <- Simulate(BaseCase)
ValidateOM(Hist_BaseCase, object = Cond_BaseCase)
Save(Hist_BaseCase, path = paste0('objects/Hist/', Name, '.hist'), overwrite = TRUE)


# ---- Robustness: [Name] ----


t = SBiomass(Hist_BaseCase)

t |> dplyr::filter(Stock=='Female') |>
  dplyr::g


library(ggplot2)
ggplot(t |> dplyr::filter( Stock=='Female'), aes(x=Year, y=Value)) +
  expand_limits(y=0) +
  geom_line()

l <- Landings(Hist_BaseCase) |>
  dplyr::filter(Sim ==2) |>
  dplyr::mutate(CalYear = floor(Year)) |>
  dplyr::group_by(CalYear) |>
  dplyr::summarise(Value = sum(Value))

ggplot(l, aes(x=CalYear, y=Value)) + geom_line() +
  expand_limits(y=0)

l |> print(n=30)
Hist_BaseCase@Landings[1,,,] |> SumOverFleet() |> SumOverStock()

Cond_BaseCase$mcmcvars[[1]]$Cmsy

2171.961+5461.863+ 6744.652+10009.384
