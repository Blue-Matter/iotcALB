library(iotcALB)

OMNames <- c('OM5a', 'OM5b', 'OM6b')

Interval    <- 3
FirstYear   <- 2000 # first assessment data year
CurrentYear <- 2020 # terminal assessment year
MPStartYear <- 2029 # first year TAC advice is set by MP

DataLag     <- 2
# 2029 TAC would use data up to 2026 or
# DataLag   = 3 # if 2026 data isn't ready in time
# (in practice calculated in 2027 and reviewed by managers in 2028
# for implementation in 2029)

mpyears <- 30 # number of years where the MP is active (ie fromMPStartYear )
pYear <- c(seq(CurrentYear + 1, by = 1, to = MPStartYear-1),
           seq(MPStartYear, by = 1, length.out = mpyears)) |> length()

# Get updated data
NewData <- ProcessNewData()

# InterimAdvice

RealCatches <- Array2DF(NewData$Landings) |>
  dplyr::filter(Year >= 2021) |>
  dplyr::rename(Mean = Value) |>
  dplyr::mutate(CV = 0, Type = 'TAC')

Assumed <- RealCatches |>
  dplyr::mutate(Year = floor(Year)) |>
  dplyr::group_by(Year, Fleet) |>
  dplyr::summarise(Mean = sum(Mean), .groups = 'drop') |>
  dplyr::group_by(Fleet) |>
  dplyr::summarise(Mean = mean(Mean),
                   .groups = 'drop') |>
  dplyr::select(Fleet, Mean) |>
  dplyr::mutate(CV = 0.1, Type = 'TAC')

FirstAssumeYear <- RealCatches$Year |> floor() |> max() + 1
AssumedYears    <- FirstAssumeYear:(MPStartYear - 1)
AssumedCatches  <- purrr::map_df(AssumedYears, \(yr)  Assumed |> dplyr::mutate(Year=yr))

InterimAdvice <- dplyr::bind_rows(RealCatches, AssumedCatches) |>
  dplyr::mutate(Stock = 'Combined')

# OM Specifications
OMSpecs <- list(
  Interval    = Interval,
  Seasons     = 4,
  FirstYear   = FirstYear,
  CurrentYear = CurrentYear,
  nYear       = length(FirstYear:CurrentYear),
  pYear       = pYear,
  MPStartYear = MPStartYear,
  DataLag     = DataLag,
  InterimAdvice = InterimAdvice,
  FleetNames  = c(paste0('LL', 1:4), 'PS', 'Other')
)



