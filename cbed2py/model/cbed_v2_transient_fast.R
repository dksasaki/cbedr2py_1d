#-----------------------------------------------------------------------------
# Transient CBED model: optimized implementation of cbed_v2_transient.R.
#
# Same equations and coefficients as cbed_v2_transient.R (including the
# generic_CBED.F90 diffusion coefficients and anammox stoichiometry); the
# right-hand sides agree to rounding (~1e-13). cbed_t_setup() is rewritten
# for speed (~19x faster transient runs), see the comments there:
#   - interface values from a precomputed spline matrix instead of
#     setup.prop.1D(xy=...) on every call
#   - tran_fast() instead of tran.1D(); valid for w >= 0 and the boundary
#     conditions used here (C.up or flux.up on top, zero gradient below)
#   - compaction formulas inline; parameters bound once
# Keep both files in step: change the equations in cbed_v2_transient.R
# first, then here, and check that the right-hand sides still agree.
#
# The rest of this header is that of cbed_v2_transient.R.
#-----------------------------------------------------------------------------
# Transient (time-evolving) CBED model with the v2 coefficients.
#
# Built from cbed_model.t() in cbed_R/cbed_v1_func.R, rewritten so that the
# reactions, rate constants and transport are those of cbed_model() in
# cbed_v2_func.R (copied next to this file).
#
# Functions
#   cbed_t_default_coefs()  v2 constants that used to be hard-coded
#   cbed_t_fill_coefs(pars) adds the default coefficients missing from pars
#   cbed_model_ss_t(pars)   steady state of the SAME equations (constant
#                           forcing from pars) -> initial condition
#   cbed_model_t(pars, times, yini, forc)  transient run
#
# Changes vs. cbed_v1_func.R's cbed_model.t():
#   - v2 rate constants: k2 = 1.5e-3, k3 = 0.9e-4, k.adj.anoxia = 0.001.
#     They and the other constants below are now parameters (defaults in
#     cbed_t_default_coefs()), so they can be changed from Python.
#   - f.OM1/f.OM2/f.OM3 come from `pars` (not hard-coded).
#   - R.TAlk uses "- 1*R.ODUox" instead of "- 2*R.ODUox".
#   - dODU scales OduDepo by svf.grid/por.grid.
#   - Compaction uses v.0 = w (as v2), not v.inf = w.
#   - Forcing comes from the `forc` data.frame argument instead of the
#     global bc.data / bc.fn.btm().
#   - The POC flux J used by k1-k3 (J.mode.k), bioturbation (J.mode.Db) and
#     bioirrigation (J.mode.irr) is chosen separately for each process:
#       0 = instantaneous forcing J.OM(t) (default)
#       1 = long-term mean, the constant pars$J.OM (as generic_CBED.F90 with
#           use_jpoc_mean_for_k / use_jpoc_mean_for_mixing = .true.)
#       2 = trailing mean of J.OM(t) over the previous J.avg.days days
#           (exact mean of the linearly interpolated forcing; before the
#           first forcing time J.OM is held at its first value)
#     The deposition flux always follows J.OM(t).
#   - switch.oxy.dependent.k is not used: the rate constants are always
#     k.adj.denit, k.adj.anoxia, k1.0 ... k3.0. v2's "switch = 0" branch
#     is k2.0 = 1.3e-4, k3.0 = 0, k.adj.denit = k.adj.anoxia = 1.
#
# Forcing (argument `forc` of cbed_model_t)
#   data.frame with a `time` column [yr] and any of
#     J.OM, O2.w, NO3.w, NH4.w, DIC.w, TAlk.w, temp, S, w
#   in CBED units. Each column is linearly interpolated in time (rule = 2:
#   held constant past the ends). Columns not given stay at their `pars`
#   value. forc = NULL -> constant forcing from pars.
#=======================================================================

library(ReacTran)
library(deSolve)
library(marelac)
library(rootSolve)

CBED_STATES  <- c("OM1","OM2","OM3","O2","NH4","NO3","ODU","DIC","TAlk")
CBED_FORCING <- c("J.OM","O2.w","NO3.w","NH4.w","DIC.w","TAlk.w","temp","S","w")


# v2 constants (defaults). Units follow the CBED code: cm, yr, mM, umol cm-2 yr-1.
cbed_t_default_coefs <- function() {
  list(
    # OM decay:  k_i = k_i.0 * J.OM^k.exp   [yr-1]
    k1.0         = 1.5e-1,
    k2.0         = 1.5e-3,
    k3.0         = 0.9e-4,
    k.exp        = 0.85,
    k.adj.denit  = 0.1,     # 1 = no effect of denitrification on decay
    k.adj.anoxia = 0.001,   # 1 = no effect of anoxia on decay
    # POC flux J used by each J.OM-dependent process:
    #   0 = instantaneous forcing J.OM(t)
    #   1 = long-term mean, the constant pars$J.OM
    #   2 = trailing mean of J.OM(t) over the previous J.avg.days days
    J.mode.k     = 0,       # k1-k3
    J.mode.Db    = 0,       # bioturbation
    J.mode.irr   = 0,       # bioirrigation
    J.avg.days   = 7,       # window of the trailing mean [d]
    # bioturbation (Archer 2002): Db = Db.coef * J.OM^Db.exp * exp(-(z/Db.xL)^2) * O2/(O2+Db.O2.half)
    Db.coef      = 0.0232,
    Db.exp       = 0.85,
    Db.xL        = 8,       # cm
    Db.O2.half   = 20,      # uM
    # bioirrigation (Archer 2002): alpha = I0 * exp(-(z/irr.xL)^2)
    irr.xL       = 1.8,     # cm
    # ODU burial (Soetaert et al. 1996): min(1, odu.bur.a * w^odu.bur.b)
    odu.bur.a    = 0.233,
    odu.bur.b    = 0.336
  )
}


# add the default coefficients that are not already in pars
cbed_t_fill_coefs <- function(pars) {
  d <- cbed_t_default_coefs()
  for (nm in names(d)) {
    if (is.null(pars[[nm]])) pars[[nm]] <- d[[nm]]
  }
  pars
}


# Trailing mean of the linear interpolant of (t, J), held constant past the
# ends (as approxfun rule = 2):   Jbar(s) = 1/tau * integral_{s-tau}^{s} J(u) du
# Exact: the integral of a piecewise-linear function is computed piecewise.
make_trailing_mean <- function(t, J, tau) {
  n  <- length(t)
  Jf <- approxfun(t, J, rule = 2)
  C  <- c(0, cumsum(diff(t) * (J[-n] + J[-1]) / 2))   # integral from t[1] to t[i]
  Cf <- function(s) {                                 # integral from t[1] to s
    s1 <- min(max(s, t[1]), t[n])
    i  <- findInterval(s1, t, rightmost.closed = TRUE)
    C[i] + (s1 - t[i]) * (J[i] + Jf(s1)) / 2 +         # inside [t[1], t[n]]
      (s - s1) * Jf(s)                                 # constant J outside
  }
  function(s) (Cf(s) - Cf(s - tau)) / tau
}


# builds the grid and the right-hand side; used by the ss and transient runs
#
# Same equations as cbed_v2_transient.R (with the F90 coefficients), written
# for speed. Everything that does not change in time is computed once here:
#   - setup.prop.1D(xy = cbind(x.mid, y)) interpolates y to the interfaces
#     with an fmm cubic spline, which is linear in y: y.int = Msp %*% y.
#     Msp ((N+1) x N) is built once by applying splinefun to unit vectors, so
#     the result equals the per-call spline up to rounding.
#   - tran.1D is replaced by tran_fast(): the same flux-form discretization
#     (C.up or flux.up at the top, zero gradient C.down = C[N] at the bottom,
#     upwind advection, AFDW = 1, A = 1) without argument checks.
#   - setup.compaction.1D is replaced by its formulas (v.0 = w).
#   - the parameters are bound once in the environment of derive(), instead
#     of with(as.list(pars)) on every call.
cbed_t_setup <- function(pars, forc = NULL) {

  # 1. forcing interpolators ------------------------------------------------
  forc_fun <- list()
  if (!is.null(forc)) {
    for (nm in intersect(CBED_FORCING, names(forc))) {
      forc_fun[[nm]] <- approxfun(forc$time, forc[[nm]], rule = 2)
    }
  }
  # value of forcing `nm` at time t (constant parameter if not forced)
  fv <- function(nm, t) {
    f <- forc_fun[[nm]]
    if (is.null(f)) pars[[nm]] else f(t)
  }

  # trailing mean of J.OM (mode 2)
  tau <- pars$J.avg.days / 365
  if (is.null(forc_fun[["J.OM"]])) {
    J.avg.fun <- function(s) pars$J.OM            # constant forcing (and steady state)
  } else if (tau > 0) {
    J.avg.fun <- make_trailing_mean(forc$time, forc$J.OM, tau)
  } else {
    J.avg.fun <- forc_fun[["J.OM"]]               # zero window = instantaneous
  }
  # POC flux used by a process, according to its mode
  J.for <- function(mode, t, J.inst) {
    switch(as.character(mode),
           "0" = J.inst,
           "1" = pars$J.OM,
           "2" = J.avg.fun(t),
           stop("J.mode.* must be 0, 1 or 2, got ", mode))
  }

  # 2. time-invariant setup (identical to cbed_v2_func.R) --------------------
  grid <- setup.grid.1D(x.up=0, x.down=pars$L, N=pars$N, dx.1=0.1)

  por.grid <- setup.prop.1D(p.exp, grid=grid, y.0=pars$por.0,
                            y.inf=pars$por.inf, x.att=pars$por.att)
  svf.grid <- setup.prop.1D(xy=cbind(grid$x.mid, 1-por.grid$mid), grid=grid)
  irr.enh  <- setup.prop.1D(func=p.exp, grid=grid, y.0=pars$irr.enh.0, y.inf=1, x.L=pars$irr.enh.L)
  phi.N    <- setup.prop.1D(xy=cbind(grid$x.mid, por.grid$mid+svf.grid$mid*pars$K.equilNH4*pars$ps), grid=grid)

  # 3. precomputed, time-invariant quantities --------------------------------
  N     <- pars$N
  dx    <- grid$dx
  dxa   <- grid$dx.aux
  xm    <- grid$x.mid
  i1 <- 1:N; i2 <- (N+1):(2*N); i3 <- (2*N+1):(3*N); i4 <- (3*N+1):(4*N)
  i5 <- (4*N+1):(5*N); i6 <- (5*N+1):(6*N); i7 <- (6*N+1):(7*N)
  i8 <- (7*N+1):(8*N); i9 <- (8*N+1):(9*N)
  # spline mid -> interfaces as a matrix (see header)
  Msp <- vapply(seq_len(N), function(j) {
    e <- numeric(N); e[j] <- 1
    splinefun(xm, e)(grid$x.int)
  }, numeric(N + 1))
  por.m  <- por.grid$mid;  por.i <- por.grid$int
  svf.m  <- svf.grid$mid;  svf.i <- svf.grid$int
  phiN.m <- phi.N$mid
  sp     <- svf.m / por.m                             # svf/por
  tortf  <- irr.enh$mid / (1 - 2*log(por.m))          # irr.enh / tortuosity
  Ktot   <- svf.m * pars$K.equilNH4 * pars$ps         # NH4 adsorption capacity
  db.shape  <- exp(-(xm / pars$Db.xL)^2)
  irr.shape <- exp(-(xm / pars$irr.xL)^2)
  por.dx <- por.m * dx
  svf.dx <- svf.m * dx

  # flux-form 1-D transport, as tran.1D (see header). VFi, VFm: volume
  # fraction at interfaces / mid; fup: imposed flux at the top (or NULL)
  tran_fast <- function(C, Cup, Di, vi, VFi, VFm, fup = NULL) {
    f <- -VFi * Di * (c(C, C[N]) - c(Cup, C)) / dxa + VFi * vi * c(Cup, C)
    if (!is.null(fup)) f[1] <- fup
    list(dC = -(f[-1] - f[-(N+1)]) / VFm / dx, flux.up = f[1], flux.down = f[N+1])
  }

  # 4. right-hand side -------------------------------------------------------
  derive <- function(times, state, parms) {

    #state variables
    OM1  <- state[i1]; OM2 <- state[i2]; OM3  <- state[i3]
    O2   <- state[i4]; NH4 <- state[i5]; NO3  <- state[i6]
    ODU  <- state[i7]; DIC <- state[i8]; TAlk <- state[i9]

    # time-varying boundary conditions
    # (these locals shadow the constant values bound from `pars`)
    J.OM   <- fv("J.OM",   times)
    O2.w   <- fv("O2.w",   times)
    NO3.w  <- fv("NO3.w",  times)
    NH4.w  <- fv("NH4.w",  times)
    DIC.w  <- fv("DIC.w",  times)
    TAlk.w <- fv("TAlk.w", times)
    temp   <- fv("temp",   times)
    w      <- fv("w",      times)

    # advection for solids (v) and solutes (u): setup.compaction.1D, v.0 = w
    v.fac <- w * (1 - por.0)
    v.m   <- v.fac / (1 - por.m); v.i <- v.fac / (1 - por.i)
    u.fac <- v.fac / (1 - por.inf) * por.inf
    u.m   <- u.fac / por.m;       u.i <- u.fac / por.i

    #define three pools of OM
    J.OM1 <- f.OM1*J.OM
    J.OM2 <- f.OM2*J.OM
    J.OM3 <- f.OM3*J.OM

    #POC flux used for reactivity, bioturbation and bioirrigation
    J.k   <- J.for(J.mode.k,   times, J.OM)
    J.Db  <- J.for(J.mode.Db,  times, J.OM)
    J.irr <- J.for(J.mode.irr, times, J.OM)

    #reactivity as function of total POC flux
    kJ <- (J.k)^k.exp
    k1 <- k1.0*kJ
    k2 <- k2.0*kJ
    k3 <- k3.0*kJ

    #bioturbation | Db as Archer 2002
    O2uM <- O2.w*1000
    Db.m <- Db.coef*(J.Db)^Db.exp * db.shape * (O2uM/(O2uM+Db.O2.half))
    Db.i <- drop(Msp %*% Db.m)

    #bioirrigation as Archer 2002, with 'if' to handle complete anoxia
    if (O2.w > 0) {
      I0 <- 11*(((atan((5*J.irr-400)/400))/pi)+0.5) - 0.9 + 20*(O2uM/(O2uM+10))*exp(-O2uM/10)*(J.irr/(J.irr+30))
    } else {
      I0 <- 0
    }
    alpha <- I0 * irr.shape

    # free-solution diffusion coefficients of generic_CBED.F90 [cm2 yr-1]
    D.O2   <- (0.031558 + 0.001428*temp)*1e4
    D.NO3  <- (0.030863 + 0.001153*temp)*1e4
    D.NH4f <- (0.030926 + 0.001225*temp)*1e4
    D.HCO3 <- (0.015179 + 0.000795*temp)*1e4
    D.H2S  <- (0.028938 + 0.001314*temp)*1e4

    # sediment diffusion (tortuosity, bioturbation) at the interfaces
    Do2.i  <- drop(Msp %*% (Db.m + D.O2  *tortf))
    Dno3.i <- drop(Msp %*% (Db.m + D.NO3 *tortf))
    Dodu.i <- drop(Msp %*% (Db.m + D.H2S *tortf))
    Ddic.i <- drop(Msp %*% (Db.m + D.HCO3*tortf))
    # ammonium (adsorption)
    Dnh4.m <- Db.m + D.NH4f*tortf
    DNH4.i <- drop(Msp %*% (por.m*Dnh4.m + Ktot*Db.m))
    wNH4.i <- drop(Msp %*% (por.m*u.m    + Ktot*v.m))

    #reactions
    # temperature dependence
    Q10.temp <- Q10^((temp-4)/10)

    fO2   <- O2/(Ks.O2+O2)
    iO2   <- Ks.O2/(Ks.O2+O2)
    fNO3  <- NO3/(Ks.NO3+NO3)
    iNO3  <- Ks.NO3/(Ks.NO3+NO3)

    #aerobic respiration
    R.O2.1 <- k1*OM1*fO2 *Q10.temp
    R.O2.2 <- k2*OM2*fO2 *Q10.temp
    R.O2.3 <- k3*OM3*fO2 *Q10.temp

    #OM oxidation by nitrate (denitrification)
    R.NO3.1 <- k.adj.denit * k1*OM1*fNO3*iO2 *Q10.temp
    R.NO3.2 <- k.adj.denit * k2*OM2*fNO3*iO2 *Q10.temp
    R.NO3.3 <- k.adj.denit * k3*OM3*fNO3*iO2 *Q10.temp

    #OM oxidation by all other electron acceptors
    R.ODU.1 <- k.adj.anoxia * k1*OM1*iO2*iNO3 *Q10.temp
    R.ODU.2 <- k.adj.anoxia * k2*OM2*iO2*iNO3 *Q10.temp
    R.ODU.3 <- k.adj.anoxia * k3*OM3*iO2*iNO3 *Q10.temp

    #total OM remineralization
    R.DIC.1 <- R.O2.1+R.NO3.1+R.ODU.1
    R.DIC.2 <- R.O2.2+R.NO3.2+R.ODU.2
    R.DIC.3 <- R.O2.3+R.NO3.3+R.ODU.3

    R.O2s  <- R.O2.1+R.O2.2+R.O2.3
    R.NO3s <- R.NO3.1+R.NO3.2+R.NO3.3
    R.ODUs <- R.ODU.1+R.ODU.2+R.ODU.3
    R.DICs <- R.DIC.1+R.DIC.2+R.DIC.3

    #nitrification
    R.NOx <- k.NOx*NH4*O2 *Q10.temp

    #anammox | NH4 + NO3 = N2
    R.ana <- k.ana*NH4*NO3 *Q10.temp

    #ODU reoxidation
    R.ODUox <- k.ODUox*ODU*O2 *Q10.temp

    #ODU burial as solid (pyrite burial) #Soetaert et al 1996
    OduDepo <- R.ODUs*min(1, odu.bur.a*(w)^odu.bur.b) *switch.odu.burial

    #TAlk generation and consumption (v2: -1*R.ODUox; F90: -0.4*R.ana)
    R.TAlk <- rNC*sp*R.O2s + (0.8+rNC)*sp*R.NO3s + (1+rNC)*sp*R.ODUs -
      2*R.NOx - 1*R.ODUox - 0.4*R.ana

    #transport
    trans.OM1  <- tran_fast(OM1,  OM1[1], Db.i,   v.i,    svf.i, svf.m, J.OM1)
    trans.OM2  <- tran_fast(OM2,  OM2[1], Db.i,   v.i,    svf.i, svf.m, J.OM2)
    trans.OM3  <- tran_fast(OM3,  OM3[1], Db.i,   v.i,    svf.i, svf.m, J.OM3)

    trans.O2   <- tran_fast(O2,   O2.w,   Do2.i,  u.i,    por.i, por.m)
    trans.NH4  <- tran_fast(NH4,  NH4.w,  DNH4.i, wNH4.i, 1,     1)
    trans.NO3  <- tran_fast(NO3,  NO3.w,  Dno3.i, u.i,    por.i, por.m)
    trans.ODU  <- tran_fast(ODU,  ODU.w,  Dodu.i, u.i,    por.i, por.m)
    trans.DIC  <- tran_fast(DIC,  DIC.w,  Ddic.i, u.i,    por.i, por.m)
    trans.TAlk <- tran_fast(TAlk, TAlk.w, Ddic.i, u.i,    por.i, por.m)

    # bioirrigation exchange
    irr.O2   <- alpha*(O2.w-O2)
    irr.NH4  <- alpha*(NH4.w-NH4)
    irr.NO3  <- alpha*(NO3.w-NO3)
    irr.ODU  <- alpha*(ODU.w-ODU)
    irr.DIC  <- alpha*(DIC.w-DIC)
    irr.TAlk <- alpha*(TAlk.w-TAlk)

    #ODEs
    dOM1 <- trans.OM1$dC - R.DIC.1
    dOM2 <- trans.OM2$dC - R.DIC.2
    dOM3 <- trans.OM3$dC - R.DIC.3

    dO2 <- trans.O2$dC - sp*R.O2s - 2*R.NOx - R.ODUox + irr.O2

    dNH4 <- 1/phiN.m*trans.NH4$dC + svf.m/phiN.m*(rNC.1*R.DIC.1 + rNC.2*R.DIC.2 + rNC.3*R.DIC.3) -
      por.m/phiN.m*(R.NOx+R.ana) + irr.NH4

    dNO3 <- trans.NO3$dC - sp*0.8*R.NO3s + R.NOx - 0.6*R.ana + irr.NO3

    dODU <- trans.ODU$dC + sp*R.ODUs - R.ODUox - sp*OduDepo + irr.ODU

    dDIC  <- trans.DIC$dC + sp*R.DICs + irr.DIC
    dTAlk <- trans.TAlk$dC + R.TAlk + irr.TAlk

    # diagnostics: the v2 set (same formulas) + a few from v1's transient
    # fluxes in umol cm-2 yr-1 (x 10/365 -> mmol m-2 d-1)
    list(c(dOM1, dOM2, dOM3, dO2, dNH4, dNO3, dODU, dDIC, dTAlk),
         o2_flux   = trans.O2$flux.up  + sum(por.dx*irr.O2),
         nh4_flux  = trans.NH4$flux.up + sum(por.dx*irr.NH4),
         no3_flux  = trans.NO3$flux.up + sum(por.dx*irr.NO3),
         tot_denit = 0.8*sum(svf.dx*R.NO3s) + 1.6*sum(por.dx*R.ana),
         anammox   = 1.6*sum(por.dx*R.ana),
         dic_flux  = trans.DIC$flux.up + sum(por.dx*irr.DIC),
         OM_flux   = trans.OM1$flux.up+trans.OM2$flux.up+trans.OM3$flux.up,
         OM_burial = trans.OM1$flux.down+trans.OM2$flux.down+trans.OM3$flux.down,
         odu_flux  = trans.ODU$flux.up + sum(por.dx*irr.ODU),
         talk_flux = trans.TAlk$flux.up + sum(por.dx*irr.TAlk),
         tot_remin = sum(svf.dx*R.DICs),
         nitrification = sum(por.dx*R.NOx),
         oxic_remin    = sum(svf.dx*R.O2s),
         OC_storage    = sum((OM1+OM2+OM3)*svf.dx)/1e6*1e4,
         Db_top    = Db.m[1],
         alpha_top = alpha[1],
         J_k   = J.k,     # POC flux used by k1-k3
         J_Db  = J.Db,    # ... by bioturbation
         J_irr = J.irr    # ... by bioirrigation
    )
  }
  # bind the (constant) parameters once: derive() sees them as variables,
  # with this setup's locals (grid, Msp, fv, ...) as the enclosing scope
  environment(derive) <- list2env(as.list(pars), parent = environment())

  list(derive = derive, grid = grid)
}


# steady state of the same equations, constant forcing taken from pars
cbed_model_ss_t <- function(pars) {
  pars <- cbed_t_fill_coefs(pars)
  s    <- cbed_t_setup(pars, forc = NULL)
  out  <- steady.1D(y=rep(0, 9*pars$N), func=s$derive, parms=pars, nspec=9,
                    names=CBED_STATES, method="runsteady")
  out$z <- s$grid$x.mid
  out
}


# transient run
cbed_model_t <- function(pars, times, yini, forc = NULL) {
  pars <- cbed_t_fill_coefs(pars)
  s    <- cbed_t_setup(pars, forc)

  out <- ode.1D(y=as.vector(yini), times=times, func=s$derive, parms=pars,
                nspec=9, names=CBED_STATES)

  # split the deSolve matrix into named pieces for Python
  N   <- pars$N
  res <- list(time = out[, 1], z = s$grid$x.mid)
  for (k in seq_along(CBED_STATES)) {
    res[[CBED_STATES[k]]] <- out[, (2 + (k-1)*N):(1 + k*N), drop = FALSE]   # (time, z)
  }
  diag_names <- colnames(out)[(2 + 9*N):ncol(out)]
  res$diag <- as.data.frame(out[, diag_names, drop = FALSE])
  res
}
