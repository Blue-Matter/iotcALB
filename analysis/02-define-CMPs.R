library(iotcALB)
library(MSEtool)


# MPs in IOTC have been generally tuned to a short period of between 11 and 15 years after the start of simulations

# Two tuning periods were considered: the standard one, from 2035 to 2039, and a longer
# period, from 2035 to 2046. Both MPs were tuned to achieve a 60% probability
# of being in the Kobe green quadrant during their respective tuning periods

# PMs:
# - how often TAC was limited by max change
# - F/FMSY
# - SB/SBMSY
# - P Green Kobe
# - SB/SB0
# - ATV
# - Prob. TAC change < 15%


# Base Case Model
OM <- readRDS('objects/OM/OM5b.om')
OM_combined <- MSEtool::CombineFleets(OM)

Hist <- readRDS('objects/Hist/OM5b.hist')
Hist_combined <- Simulate(OM_combined)

PlotSBiomass(Hist, Stocks = 'Female', Season = 1)
PlotSBiomass(Hist_combined, Stocks = 'Female', Season = 1)

PlotLandings(Hist, byStock = 'sum', byFleet = FALSE, AggregateYear = TRUE)
PlotLandings(Hist_combined, byStock = 'sum', byFleet = FALSE, AggregateYear = TRUE)


MSE <- Project(Hist_combined,
               MPs = c(SP1 = MSEtool::SetMPArgs(SurplusProduction,
                                                EstDepletion = TRUE),
                       IR1 = IndexRate,
                       IT1 = IndexTarget),
               nSim = 50
)

MSE_test <- Project(Hist_combined, MPs = c('NoFishing', 'CurrentCatch'))

plot(MSE_test)

plot(MSE)

PlotSBiomass(MSE, Season = 1, Stock = 'Female', byMP = TRUE)
PlotSBiomass(MSE, Season = 1, relative = 'BMSY', Stock = 'Female', byMP = TRUE)


PlotLandings(MSE, byStock = 'sum', byFleet = FALSE, byMP = TRUE, AggregateYear = TRUE)

sims <- 1
plot(apply(MSE@SBiomass[sims,1,,1, drop=FALSE], 3, mean), type='l')


Hist_combined@OM@StockTargeting@Targeting |> dim()

Log(MSE)

df = SB_SBMSY(MSE) |> dplyr::filter(Stock =='Female', Sim ==14)

r = SB_SB0(MSE) |> dplyr::filter(Sim == 14, Stock == 'Female') |>
  dplyr::mutate(Year = floor(Year)) |>
  dplyr::group_by(Year) |>
  dplyr::summarise(Value = mean(Value))

plot(r$Year, r$Value, type='l', ylim=c(0,1))

OM_combined@InterimAdvice |>
  dplyr::mutate(Year = floor(Year)) |>
  dplyr::group_by(Year) |>
  dplyr::summarise(Value = sum(Mean))

Hist_combined |> MSYRefs() |> MSYLandings() |> Array2DF() |> dplyr::filter(Sim == 14)
Hist_combined |> MSYRefs() |> MSYLandings() |> hist(breaks = 20)
Hist_combined |> MSYRefs() |> MSYLandings() |> quantile(c(0.05, 0.5, 0.95))

sb0 <- SB0(Hist) |> Array2DF() |> dplyr::filter(Stock == 'Female', Year == 2000)
sb0 |> dplyr::filter(Sim == 14)
sb0 |> dplyr::arrange(Value) |> head(10)

hist(sb0$Value, breaks = 100)


t = df |> dplyr::filter(Sim == 14)
plot(t$Year, t$Value, type='l')

n <- Number(MSE) |> dplyr::filter(Sim == 14)
plot(n$Year, n$Value)

d <- c(MSE@OM@Stock$Female@SRR@RecDevHist[14,],MSE@OM@Stock$Female@SRR@RecDevProj[14,])
d <- d[d>0]
plot(d, type='b')

l <- Removals(MSE)
ldf <- l |> dplyr::filter(Sim == 14, MP == 'IR1') |> dplyr::group_by(Year) |>
  dplyr::summarise(Value = sum(Value))


plot(ldf$Year, ldf$Value, type='l')



# Tuning targets:
# - maximise mean yield from years 11 - end of projections
# - maintain at least 60% probability of green kobe in years 2035 - 2046



# Tune with 80 simulations

# Full run with all 500

Hist_combined@OM@Obs$Combined$LL1@Survey


Data <- DataList$`1`$Combined

# ---- Index Rate MP ----

# ---- Index Target MP ----

# ---- Surplus Production Model ----
Data |> CPUE()

dnames <- dimnames(Data@Landings@Value)

Data@Years <- as.numeric(dnames$Year)

fl <- tempfile()
fl
saveRDS(Data, fl)
t <- SurplusProduction(Data,
                       IndexSource = 'CPUE',
                       Diagnostics = 'full',
                       EstDepletion = TRUE)

t@Misc$SurplusProduction$Summary




t@Misc$SurplusProduction$par


