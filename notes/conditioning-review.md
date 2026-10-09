# Review of the ABC Conditioning (albMSE): Detailed Notes

Reference notes for the review of the Indian Ocean albacore conditioning model
(Hillary & Mosqueira), compiled while importing the OMs into `iotcALB`
(October 2026). The TSD (`tsd/operating-models.qmd`, *Review of the
Conditioning*) has a short summary; these are the details.

- Source: local clone of <https://github.com/iagomosqueira/albMSE>, branch
  `main`, commit `e29a97f` (2025-12-05), up to date with GitHub on 2026-10-07.
- Evidence comes from the scripts (`data_abc5a.R`, `data_abc5b.R`,
  `data_abc6b.R`, `utilities.R`) and from the ABC output tracked in the repo
  (`data/omX/image_abcX.rda`: the sampler settings saved before the run;
  `data/omX/mcmc_abcX.rda`: the MCMC draws; `data/omX/mcvars_abcX.rda`: the
  posterior quantities imported into `iotcALB` as `CondData_OMX`).
- **Caveat**: the repository may not contain the exact code or outputs used for
  the results in the IOTC documents. Everything below should be confirmed with
  the authors before it is presented as a problem with their work.
- MSY check: `analysis/99-check-abc-msy.R` (needs the albMSE clone at
  `../albMSE`).

---

## 1. Specification of the robustness OMs

Described specifications (Mosqueira & Hillary 2025, IOTC-2026-WPM17(MSE)-03):

| OM   | Described |
|------|-----------|
| OM5b | Base Case: LL1 (NW) CPUE |
| OM5a | Robustness: SW CPUE (called LL2 (SW) in WPM17; LL3 (SW) in TCMP10-04 and in the 2022 assessment fleet structure) |
| OM6b | Robustness: LL1 (NW) CPUE with a 1% annual increase in LL catchability |

### OM5a: conditioned on LL1, not SW

- `data_abc5a.R` line 126: `fcpue <- 3` (the third CPUE series, SW). Line 132:
  `fcpue <- 1`, which overwrites it. The first block (lines 126--130) does not
  calculate `sd.cpue` (no loop), so it looks like a leftover.
- Apart from file names, the only differences from `data_abc5b.R` that affect
  the fit are that line and `rwsd[paridx[[3]]]` (0.04 vs 0.025, the random-walk
  SD of the selectivity parameters).
- The saved image `data/om5a/image_abc5a.rda` has `fcpue = 1`.
- The SW index appears only in the FLR observation model for the projections
  (line 295, `base$ids[9:12]`), which `iotcALB` does not use.

### OM6b: catchability trend set but not used

- `data_abc6b.R` lines 153--155: `qtrend <- TRUE`, `delq <- 0.01`,
  `qt <- exp(delq*(0:(ny-1))) %o% rep(1,4)`. The saved image has
  `qtrend = TRUE`, `delq = 0.01`.
- All three OMs are fitted with `mcmc5.abc()` (`utilities.R` lines
  1612--1853), which never uses `qtrend` or `qt`; neither does `sim()`. The
  CPUE discrepancy is `log(I[,,fcpue]/xx$I)` with no trend.
- The trend is applied (as `log(I[,,fcpue]/(xx$I*qt))`) only in the older
  samplers `mcmc2.abc()`, `mcmc2a.abc()`, `mcmc3.abc()`, `mcmc3a.abc()`,
  `mcmc4.abc()` (e.g. line 457) and in `plot.mcmc.cpue()`.
- The observation-model code of `data_abc6b.R` is identical to OM5b, so no
  trend in the FLR projections either.

### Empirical check on the ABC output

Median predicted annual log-CPUE vs the observed annual log-CPUE of each series
(`alb_abcdata.rda`, `I[year, season, series]`):

| OM   | cor. series 1 (LL1, NW) | cor. series 3 (SW) | RMSE vs 1 | RMSE vs 3 |
|------|------|------|------|------|
| OM5b | 0.90 | 0.48 | 0.10 | 0.18 |
| OM5a | 0.89 | 0.47 | 0.10 | 0.18 |
| OM6b | 0.89 | 0.48 | 0.10 | 0.18 |

A model conditioned on SW would track series 3. A 1%/yr catchability increase
over 2000--2020 would make the fitted biomass trend about 0.2 log units flatter
than the CPUE. Neither is seen.

Posterior summaries (median, 5th--95th percentiles):

| OM   | R0 (millions)     | SB/SB0 2000      | SB/SB0 2010      | SB/SB0 2020      |
|------|-------------------|------------------|------------------|------------------|
| OM5b | 14.7 (9.4--26.5)  | 0.76 (0.64--0.85) | 0.61 (0.46--0.71) | 0.42 (0.25--0.53) |
| OM5a | 13.5 (10.8--19.6) | 0.74 (0.58--0.86) | 0.59 (0.44--0.72) | 0.42 (0.27--0.53) |
| OM6b | 13.7 (10.3--20.3) | 0.71 (0.59--0.83) | 0.58 (0.42--0.71) | 0.42 (0.24--0.54) |

Median h (0.80) and M (0.075) are the same in all three. Median SSB is about
7--10% higher in OM5b (88 vs 81--82 kt in 2000), from the upper tail of R0
(absolute size is the least identified quantity; the CPUE only informs the
trend). The differences are consistent with separate, short MCMC chains (500
retained draws from 10 parallel chains) with the same likelihood.

**Conclusion**: in the repository code and output, OM5a and OM6b are
replicates of OM5b. Only OM5b is projected in the current MSE
(`analysis/03-run-projections.R`, `ProjectOMs`).

---

## 2. MSY calculation

### The code (`utilities.R`)

1. `msyfn()` (line 250) evaluates
   `msypdyn(c(ns,na,nf), srec, R0, h, psi, M, ...)`. `R0`, `h`, `M` are not
   arguments. R resolves free variables in the environment where the function
   was **defined** (lexical scoping), not where it is called. `msyfn()` is
   defined at the top level, so it uses the global objects: `R0 = 14e6`,
   `h = 0.8` (`data_abc5b.R` lines 44 and 46, the sampler's starting values) and
   `M = 0.075` (`boot/data/alb_abcdata.rda`).
2. `sim()` receives the draw's `R0`, `h`, `M` as arguments and calls
   `optimise(msyfn, interval = c(0, 0.9), ...)` (line 1894). The arguments of
   `sim()` are not visible to `msyfn()`, so `Hmsy` and `Cmsy` are those of the
   global values. `Bmsy` is then computed with `msypdyn()` using the draw's
   parameters (line 1897) but at the global-value `Hmsy`.
3. `mcmc5.abc()` uses `Bmsy` and `Hmsy` in the status priors (lines 1650--1692
   and 1756--1799): SB/SBMSY in the last two years (`ybmsy = c(ny-1, ny)`,
   means 2.25 and 2, SD 0.35) and an overfishing penalty on HR/HRMSY in every
   year (`yof = 1:ny`, SD 0.5). Both enter the acceptance step. `Cmsy` is only
   reported. The depletion prior (SB/B0 in year 1, mean 0.5, SD 0.1) uses B0
   and is not affected.
4. `objfn.init()` (line 235), used by `sim()` to solve the initial harvest
   rates, also uses the global `M`.

Minimal illustration:

```r
f <- function() x
g <- function(x) f()
x <- 1
g(5)  # returns 1, not 5
```

### Re-running `sim()` (150 posterior draws of OM5b)

| Quantity (5th; 50th; 95th) | As written | `msyfn()` with the draw's R0, h, M |
|---|---|---|
| MSY (t) | 46,600; 47,400; 48,100 | 31,500; 52,700; 98,000 |
| cor(MSY, R0) | -0.26 | 0.97 |
| Harvest rate at MSY | 0.46; 0.61; 0.90 | 0.45; 0.62; 0.90 |
| SB/SBMSY, last year | 1.22; 2.02; 2.49 | 1.14; 1.98; 2.81 |

- Posterior h is 0.73--0.86 and M 0.069--0.081, close to the globals, so the
  effect on `Hmsy` is small (change correlated 0.95 with h).
- SB/SBMSY changes per draw (correlation 0.82 between the two calculations).
- Also fixing `objfn.init()`: initial depletion changes by < 0.01, but the CPUE
  log-likelihood changes by up to 17 units for one draw.

### Effect on the conditioning (approximate)

- Status prior log-density changes by -5.3 to +1.5 (5th--95th) per draw.
- Importance reweighting of the accepted draws by the corrected target density:
  effective sample size 53 of 102.
- Reweighted posterior: SB/SBMSY in the last year narrower (5th percentile 1.76
  vs 1.55; 95th 2.53 vs 2.86); MSY similar (median 59 vs 57 kt).
- Weights uncorrelated with MSY (0.07): the fix would trim the low-status tail
  but not remove the low-productivity simulations.
- Only re-running the conditioning gives the real effect.

### Consequences in the current MSE

- ABC `Cmsy` is ~constant (46.6--48.2 kt); OM MSY (MSEtool `CalcMSY()`) ranges
  32--91 kt and is proportional to B0 (correlation 0.99).
- With the OM MSY, 2021--2023 catches exceed MSY in about 22% of simulations;
  the 15 of 500 simulations that collapse before 2029 have MSY 26--33 kt and
  2021--2023 catch/MSY of 1.3--1.6.
- The MSE uses the OM reference points throughout.

### Suggested correction

```r
msyfn <- function(H, ph, sela, R0, h, M) { ... }
msy <- optimise(msyfn, interval = c(0, 0.9), ph = ph, sela = sela,
                R0 = R0, h = h, M = M, maximum = TRUE)
```

and the same for `objfn.init()` (pass `M`). Check the upper bound of the
search interval (below).

---

## 3. Other observations (to confirm with the authors)

- **Stored draws that fail the length composition criterion.** 48 of 150
  stored draws have KL(LF) > `KLmax` (0.8) when re-evaluated with their stored
  h and M. In `mcmc5.abc()`, (h, M) is redrawn from its prior with probability
  `acphmu` = 0.25 each iteration (line 1706, `rmvnorm(1, c(hmu, Mmu), Sigma)`)
  with no accept/reject step. The stored draw then combines parameters accepted
  under the previous (h, M) with the new (h, M), and the next acceptance step
  compares the discrepancy under the new (h, M) with `dtotold` computed under
  the previous (h, M). This may be intentional (h and M fixed to their prior,
  a "cut"), but it means some stored parameter sets were never evaluated
  together.
- **Hmsy at the bound.** In about 9% of draws the harvest rate at MSY is at
  the upper bound of the `optimise()` interval (0.9), so MSY is not found.
- **Iteration 13 of OM5b.** `data_abc5b.R` replaces `mcvars[[13]]` with
  `mcvars[[113]]` (comment `# BUG iter 13`), so OM5b has a duplicated draw.

---

## 4. Data used for conditioning vs the 2025 assessment data

- Conditioning data: 2022 assessment, 2000--2020. Updated data: 2025
  assessment, 2000--2023, with revised historical catches and a revised CPUE
  standardisation.
- **CPUE**: the updated LL1 index declines more steeply over the last decade;
  2021--2023 values average 34% below the 2016--2020 mean. The OMs do not
  predict a decline this large in most simulations.
- **Catch** (annual totals, all fleets):
  - mean 2000--2020: 36,600 t (conditioning) vs 38,200 t (updated), +4%;
  - most years within 10%; largest revisions 2013 (+52%, 32,700 vs 49,800 t),
    2014 (+15%) and 2015 (+11%); 2020 -10%;
  - 2016--2020 means almost the same (39,400 vs 39,200 t);
  - by fleet: "Other" fleet revised up throughout 2000--2016 (mean 1,700 vs
    4,100 t); catches reallocated among LL fleets (LL2 higher in 2011--2018,
    lower in 2019--2020; LL4 lower in 2014--2017; LL1 lower in 2000--2001 and
    2018--2020);
  - 2021--2023 updated totals: 34,500, 48,400, 41,800 t (the 2022 peak mostly
    LL4).

---

## 5. Harvest rate vs Baranov

The ABC model uses pulse fishing (harvest rates); the OMs use the Baranov
equation, with F = -log(1 - H). The OMs track the ABC model within about 5% in
numbers for most simulations, up to ~20% for a few high-exploitation
simulations by 2020. See the TSD (*Differences Between the ABC Conditioning
Model and the OMs*).

---

## 6. Recommendations

1. Recondition on the 2025 assessment data (catch and CPUE to 2023, revised
   CPUE standardisation), which are the data the MPs will use.
2. Correct the MSY calculation (and `objfn.init()`); confirm the (h, M)
   sampling; widen or check the Hmsy search interval.
3. Implement the robustness OMs as intended: SW CPUE in the likelihood (OM5a);
   a sampler that applies `qtrend` (OM6b); and the climate change and
   recruitment OMs (OM5b_cc, OM5b_rec), which are not in the repository.
4. Consider conditioning with the Baranov equation so the OMs reproduce the
   conditioning model exactly.
5. Re-tune the CMPs on the reconditioned OMs.
