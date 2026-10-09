# Check of the MSY calculation in the ABC conditioning model (albMSE, OM5b)
#
# albMSE utilities.R msyfn() evaluates the equilibrium catch with R0, h, and M
# found in the global environment (R0 = 14e6 and h = 0.8 in data_abc5b.R, M =
# 0.075 from boot/data/alb_abcdata.rda), not the values of the draw passed to
# sim(). This script re-runs sim() for posterior draws of OM5b:
#   A: code as written
#   B: msyfn() uses the draw's R0, h, and M
#   C: B, and objfn.init() (initial harvest rates) uses the draw's M
# and compares the MSY quantities, the stock status priors of mcmc5.abc(), and
# the change in the ABC target density (importance weights).
#
# Requires a local clone of https://github.com/iagomosqueira/albMSE with the
# OM5b outputs (data/om5b), Rcpp, and mvtnorm. Not part of the MSE workflow.
AlbMSEDir <- '../albMSE'
nDraws    <- 150

owd <- setwd(AlbMSEDir)
suppressMessages({library(Rcpp); library(mvtnorm)})
load('data/om5b/image_abc5b.rda')
load('data/om5b/mcmc_abc5b.rda')
sourceCpp('utilities/init_pdyn.cpp'); sourceCpp('utilities/msy_pdyn.cpp'); sourceCpp('utilities/pdyn_lfcpue.cpp')
setwd(owd)
cat('Global values used by msyfn(): R0 =', R0, ' h =', h, ' M =', M, '\n')
R0g <- R0; hg <- h; Mg <- M
msyfn_asis <- msyfn
msyfn_draw <- function(H, ph, sela) {
  resx <- msypdyn(c(ns,na,nf), srec, R0f, hf, psi, Mf, as.vector(mata), as.vector(wta), as.vector(sela), H * ph)
  sum(resx$C)
}

Pars <- do.call(rbind, lapply(mczzz, `[[`, 'pars'))
set.seed(1); Draws <- sample(nrow(Pars), nDraws)

# discrepancy terms of mcmc5.abc()
Disc <- function(xx, epsr, sigR) {
  phat <- xx$LF; kllf <- pobs * log(pobs / phat)
  dlf <- sum(apply(kllf, 2, function(x) sum(x[!is.nan(x)])))
  resq <- log(I[,,fcpue] / xx$I); lnq <- apply(resq, 2, mean)
  resq <- t(apply(resq, 1, function(x) x - lnq))
  dcpue <- sum(dnorm(resq, 0, sdcpue, TRUE))
  Bmsyrat <- (xx$S[,3] / xx$Bmsy)[ybmsy]
  dSSB <- (xx$S[,srec-1] / xx$B0)[ydep]
  hy <- apply(xx$H[yof,,], c(1,2), sum)
  hmsyrat <- apply(apply(hy, 1, function(x) x / xx$Hmsy), 2, mean)
  zof <- pmax(hmsyrat - 1, 0)
  pB <- sum(dnorm(Bmsyrat, mubmsy, sdbmsy, TRUE)); pD <- dnorm(dSSB, mudep, sddep, TRUE)
  pO <- sum(dnorm(zof, 0, sdof, TRUE) - dnorm(0, 0, sdof, TRUE))
  c(Cmsy = xx$Cmsy, Hmsy = sum(xx$Hmsy), Bmsy_B0 = xx$Bmsy / xx$B0, Bmsyrat_ny = Bmsyrat[2],
    hmsyrat_max = max(hmsyrat), dep1 = dSSB, pB = pB, pD = pD, pO = pO, dcpue = dcpue, dlf = dlf,
    dtot = dcpue + pB + pD + pO + sum(dnorm(epsr, 0, sigR, TRUE)) - dlf)
}

Out <- lapply(Draws, function(i) {
  p <- Pars[i, ]
  R0x <- exp(p[1]); depx <- ilogit(p[2]); epsrx <- p[3:(ny+1)]; selvx <- exp(p[(ny+2):npar])
  selparsx <- cbind(selvx[1:nselg], selvx[(nselg+1):(2*nselg)], selvx[(2*nselg+1):(3*nselg)])
  hx <- p[npar+1]; Mx <- p[npar+2]; sigR <- p[npar+3]
  R0f <<- R0x; hf <<- hx; Mf <<- Mx
  run <- function(variant) {
    R0 <<- R0g; h <<- hg; M <<- Mg; msyfn <<- msyfn_asis
    if (variant %in% c('B', 'C')) msyfn <<- msyfn_draw
    if (variant == 'C') { R0 <<- R0x; h <<- hx; M <<- Mx }
    Disc(sim(R0x, depx, hx, Mx, selparsx, epsrx, dms, pctarg, selidx), epsrx, sigR)
  }
  rbind(A = run('A'), B = run('B'), C = run('C'))
})
R0 <- R0g; h <- hg; M <- Mg; msyfn <- msyfn_asis
A <- t(sapply(Out, `[`, 'A', )); B <- t(sapply(Out, `[`, 'B', )); C <- t(sapply(Out, `[`, 'C', ))
hM <- Pars[Draws, npar + 1:2]; colnames(hM) <- c('h', 'M')
R0d <- exp(Pars[Draws, 1])

q <- function(x, p = c(.05, .5, .95)) round(quantile(x, p), 3)
cat('\nh and M of the draws (5/50/95%):\n'); print(apply(hM, 2, q))
cat('\nQuantiles (5/50/95%), A: as written | B: msyfn fixed | C: msyfn and objfn.init fixed\n')
for (v in c('Cmsy', 'Hmsy', 'Bmsy_B0', 'Bmsyrat_ny', 'hmsyrat_max', 'dep1', 'dlf', 'dcpue'))
  cat(sprintf('%-12s', v), 'A', q(A[,v]), '| B', q(B[,v]), '| C', q(C[,v]), '\n')
cat('\ncor(Cmsy, R0): A', round(cor(A[,'Cmsy'], R0d), 2), ' B', round(cor(B[,'Cmsy'], R0d), 2), '\n')
cat('cor(Hmsy B/A, h):', round(cor(B[,'Hmsy'] / A[,'Hmsy'], hM[,1]), 2), '\n')
cat('Hmsy at the upper bound of optimise() (0.9): A', mean(A[,'Hmsy'] > 0.899), '\n')

# Draws that pass the length composition criterion as stored
ok <- A[,'dlf'] < KLmax
cat('\nStored draws with KL(LF) < KLmax when re-evaluated:', sum(ok), 'of', length(ok), '\n')

# Importance weights of the corrected target density (B) relative to the as-written one (A)
lw <- (B[,'dtot'] - A[,'dtot'])[ok]; w <- exp(lw - max(lw)); w <- w / sum(w)
cat('Change in the status prior (log), B - A:', q((B[,'pB'] + B[,'pO']) - (A[,'pB'] + A[,'pO'])), '\n')
cat('Effective sample size of the reweighted draws:', round(1 / sum(w^2)), 'of', sum(ok), '\n')
wq <- function(x, w, p) { o <- order(x); approx(cumsum(w[o]), x[o], p, rule = 2)$y }
P <- c(.05, .25, .5, .75, .95)
for (v in c('Cmsy', 'Bmsyrat_ny', 'hmsyrat_max'))
  cat(sprintf('%-12s', v), 'unweighted', round(quantile(B[ok, v], P), 2), '| reweighted', round(wq(B[ok, v], w, P), 2), '\n')
cat('cor(log weight, corrected Cmsy):', round(cor(lw, B[ok, 'Cmsy']), 2), '\n')
