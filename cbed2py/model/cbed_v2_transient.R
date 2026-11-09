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
#   - k1-k3 and Db/bioirrigation use the instantaneous J.OM(t) by default.
#     k.use.J.mean = 1 / mix.use.J.mean = 1 use the long-term mean J.OM
#     (the constant pars$J.OM) instead, as generic_CBED.F90 does with
#     use_jpoc_mean_for_k / use_jpoc_mean_for_mixing = .true. (its defaults).
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
    k.use.J.mean = 0,       # 1: k from the constant pars$J.OM (long-term mean, as
                            #    generic_CBED.F90 with use_jpoc_mean_for_k = .true.)
                            # 0: k from the instantaneous forcing J.OM(t)
    mix.use.J.mean = 0,     # 1: Db and bioirrigation from the constant pars$J.OM
                            #    (generic_CBED.F90, use_jpoc_mean_for_mixing = .true.)
                            # 0: from the instantaneous forcing J.OM(t)
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


# builds the grid and the right-hand side; used by the ss and transient runs
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
    if (is.null(forc_fun[[nm]])) pars[[nm]] else forc_fun[[nm]](t)
  }

  # 2. time-invariant setup (identical to cbed_v2_func.R) --------------------
  grid <- setup.grid.1D(x.up=0, x.down=pars$L, N=pars$N, dx.1=0.1)

  por.grid <- setup.prop.1D(p.exp, grid=grid, y.0=pars$por.0,
                            y.inf=pars$por.inf, x.att=pars$por.att)
  svf.grid <- setup.prop.1D(xy=cbind(grid$x.mid, 1-por.grid$mid), grid=grid)
  irr.enh  <- setup.prop.1D(func=p.exp, grid=grid, y.0=pars$irr.enh.0, y.inf=1, x.L=pars$irr.enh.L)
  phi.N    <- setup.prop.1D(xy=cbind(grid$x.mid, por.grid$mid+svf.grid$mid*pars$K.equilNH4*pars$ps), grid=grid)

  # 3. right-hand side -------------------------------------------------------
  derive <- function(times, state, pars) {
    with(as.list(pars), {

      #state varibles
      OM1      <- state[(0*N+1):(1*N)]
      OM2      <- state[(1*N+1):(2*N)]
      OM3      <- state[(2*N+1):(3*N)]
      O2       <- state[(3*N+1):(4*N)]
      NH4      <- state[(4*N+1):(5*N)]
      NO3      <- state[(5*N+1):(6*N)]
      ODU      <- state[(6*N+1):(7*N)]  # ODU = oxygen demand unit
      DIC      <- state[(7*N+1):(8*N)]
      TAlk     <- state[(8*N+1):(9*N)]

      J.OM.mean <- J.OM      # constant value from pars (long-term mean)

      # time-varying boundary conditions
      # (these locals shadow the constant values coming from `pars`)
      J.OM   <- fv("J.OM",   times)
      O2.w   <- fv("O2.w",   times)
      NO3.w  <- fv("NO3.w",  times)
      NH4.w  <- fv("NH4.w",  times)
      DIC.w  <- fv("DIC.w",  times)
      TAlk.w <- fv("TAlk.w", times)
      temp   <- fv("temp",   times)
      S      <- fv("S",      times)
      w      <- fv("w",      times)

      # advection for solids (v) and solutes (u), as v2 (v.0 = w)
      comp   <- setup.compaction.1D(v.0=w, por.0=por.0, por.inf=por.inf,
                                    por.grid=por.grid)
      u.grid <- comp$u
      v.grid <- comp$v

      #define three pools of OM
      J.OM1 <- f.OM1*J.OM
      J.OM2 <- f.OM2*J.OM
      J.OM3 <- f.OM3*J.OM

      #reactivity as function of total POC flux
      J.k <- if (k.use.J.mean == 1) J.OM.mean else J.OM
      k1 <- k1.0*(J.k)^k.exp
      k2 <- k2.0*(J.k)^k.exp
      k3 <- k3.0*(J.k)^k.exp

      #POC flux used for bioturbation and bioirrigation
      J.mix <- if (mix.use.J.mean == 1) J.OM.mean else J.OM

      #bioturbation | Db as Archer 2002
      Db.0 <- Db.coef*(J.mix)^Db.exp
      Db.grid <- setup.prop.1D(xy = cbind(x=grid$x.mid,
                                          y=(Db.0*exp(-(grid$x.mid/Db.xL)^2)*((O2.w*1000)/(O2.w*1000+Db.O2.half)))
      ), grid=grid)

      #bioirrigation as Archer 2002, with 'if' to handle complete anoxia
      if (O2.w > 0) {
        I0 <- 11*(((atan((5*J.mix-400)/400))/pi)+0.5) - 0.9 + 20*((O2.w*1000)/(O2.w*1000+10))*exp(-O2.w*1000/10)*(J.mix/(J.mix+30))
      } else {
        I0 <- 0
      }
      alpha <- setup.prop.1D(xy=cbind(grid$x.mid,
                                      I0*exp(-(grid$x.mid/irr.xL)^2)),
                             grid = grid)

      # diffusion terms
      D <- diffcoeff(S=S, t=temp, species=c("O2","NO3","NH4","Fe","SO4","H2S","CH4","Mn","HCO3"))*100^2*365*86400

      # Correct for tortuosity
      Do2.grid  <- setup.prop.1D(xy=cbind(grid$x.mid,
                                          Db.grid$mid+(D$O2*irr.enh$mid)/(1-2*log(por.grid$mid))),
                                 grid=grid)
      Dno3.grid <- setup.prop.1D(xy=cbind(grid$x.mid,
                                          Db.grid$mid+(D$NO3*irr.enh$mid)/(1-2*log(por.grid$mid))),
                                 grid=grid)
      Dnh4.grid <- setup.prop.1D(xy=cbind(grid$x.mid,
                                          Db.grid$mid+(D$NH4*irr.enh$mid)/(1-2*log(por.grid$mid))),
                                 grid=grid)
      Dodu.grid <- setup.prop.1D(xy=cbind(grid$x.mid,
                                          Db.grid$mid+(D$H2S*irr.enh$mid)/(1-2*log(por.grid$mid))),
                                 grid=grid)
      Ddic.grid <- setup.prop.1D(xy=cbind(grid$x.mid,
                                          Db.grid$mid+(D$HCO3*irr.enh$mid)/(1-2*log(por.grid$mid))),
                                 grid=grid)
      # ammonium (adsorption)
      D.NH4 <- setup.prop.1D(xy=cbind(grid$x.mid, por.grid$mid*Dnh4.grid$mid + svf.grid$mid*K.equilNH4*ps*Db.grid$mid),
                             grid=grid)
      w.NH4 <- setup.prop.1D(xy=cbind(grid$x.mid, por.grid$mid*u.grid$mid + svf.grid$mid*K.equilNH4*ps*v.grid$mid),
                             grid=grid)

      #reactions
      # temperature dependence
      Q10.temp <- Q10^((temp-4)/10)

      #aerobic respiration
      R.O2.1 <- k1*OM1*(O2/(Ks.O2+O2)) *Q10.temp
      R.O2.2 <- k2*OM2*(O2/(Ks.O2+O2)) *Q10.temp
      R.O2.3 <- k3*OM3*(O2/(Ks.O2+O2)) *Q10.temp

      #OM oxidation by nitrate (denitrification)
      R.NO3.1 <- k.adj.denit * k1*OM1*(NO3/(Ks.NO3+NO3))*(Ks.O2/(Ks.O2+O2)) *Q10.temp
      R.NO3.2 <- k.adj.denit * k2*OM2*(NO3/(Ks.NO3+NO3))*(Ks.O2/(Ks.O2+O2)) *Q10.temp
      R.NO3.3 <- k.adj.denit * k3*OM3*(NO3/(Ks.NO3+NO3))*(Ks.O2/(Ks.O2+O2)) *Q10.temp

      #OM oxidation by all other electron acceptors
      R.ODU.1 <- k.adj.anoxia * k1*OM1*(Ks.O2/(Ks.O2+O2))*(Ks.NO3/(Ks.NO3+NO3)) *Q10.temp
      R.ODU.2 <- k.adj.anoxia * k2*OM2*(Ks.O2/(Ks.O2+O2))*(Ks.NO3/(Ks.NO3+NO3)) *Q10.temp
      R.ODU.3 <- k.adj.anoxia * k3*OM3*(Ks.O2/(Ks.O2+O2))*(Ks.NO3/(Ks.NO3+NO3)) *Q10.temp

      #total OM remineralization
      R.DIC.1 <- R.O2.1+R.NO3.1+R.ODU.1
      R.DIC.2 <- R.O2.2+R.NO3.2+R.ODU.2
      R.DIC.3 <- R.O2.3+R.NO3.3+R.ODU.3

      #nitrification
      R.NOx <- k.NOx*NH4*O2 *Q10.temp

      #anammox | NH4 + NO3 = N2
      R.ana <- k.ana*NH4*NO3 *Q10.temp

      #ODU reoxidation
      R.ODUox <- k.ODUox*ODU*O2 *Q10.temp

      #ODU burial as solid (pyrite burial) #Soetaert et al 1996
      OduDepo <- (R.ODU.1+R.ODU.2+R.ODU.3)*min(1, odu.bur.a*(w)^odu.bur.b) *switch.odu.burial

      #TAlk generation and consumption (v2: -1*R.ODUox)
      R.TAlk <- rNC*svf.grid$mid/por.grid$mid*(R.O2.1+R.O2.2+R.O2.3) +
        (0.8+rNC)*svf.grid$mid/por.grid$mid*(R.NO3.1+R.NO3.2+R.NO3.3) +
        (1+rNC)*svf.grid$mid/por.grid$mid*(R.ODU.1+R.ODU.2+R.ODU.3) -
        2*R.NOx - 1*R.ODUox

      #transport
      trans.OM1 <- tran.1D(C=OM1, flux.up = J.OM1, D=Db.grid, v=v.grid, dx=grid,VF=svf.grid)
      trans.OM2 <- tran.1D(C=OM2, flux.up = J.OM2, D=Db.grid, v=v.grid, dx=grid,VF=svf.grid)
      trans.OM3 <- tran.1D(C=OM3, flux.up = J.OM3, D=Db.grid, v=v.grid, dx=grid,VF=svf.grid)

      trans.O2   <- tran.1D(C=O2,   C.up=O2.w,   D=Do2.grid,  v=u.grid, dx = grid, VF=por.grid)
      trans.NH4  <- tran.1D(C=NH4,  C.up=NH4.w,  D=D.NH4,     v=w.NH4,  dx = grid)
      trans.NO3  <- tran.1D(C=NO3,  C.up=NO3.w,  D=Dno3.grid, v=u.grid, dx = grid, VF=por.grid)
      trans.ODU  <- tran.1D(C=ODU,  C.up=ODU.w,  D=Dodu.grid, v=u.grid, dx = grid, VF=por.grid)
      trans.DIC  <- tran.1D(C=DIC,  C.up=DIC.w,  D=Ddic.grid, v=u.grid, dx = grid, VF=por.grid)
      trans.TAlk <- tran.1D(C=TAlk, C.up=TAlk.w, D=Ddic.grid, v=u.grid, dx = grid, VF=por.grid)

      #ODEs
      dOM1 <- trans.OM1$dC - R.O2.1 - R.NO3.1 - R.ODU.1
      dOM2 <- trans.OM2$dC - R.O2.2 - R.NO3.2 - R.ODU.2
      dOM3 <- trans.OM3$dC - R.O2.3 - R.NO3.3 - R.ODU.3

      dO2 <- trans.O2$dC - svf.grid$mid/por.grid$mid*(R.O2.1+R.O2.2+R.O2.3) - 2*R.NOx - R.ODUox + alpha$mid*(O2.w-O2)

      dNH4 <- 1/phi.N$mid*trans.NH4$dC + svf.grid$mid/phi.N$mid*(rNC.1*R.O2.1+rNC.2*R.O2.2 + rNC.3*R.O2.3 +
                                                                   rNC.1*R.NO3.1+rNC.2*R.NO3.2+rNC.3*R.NO3.3+
                                                                   rNC.1*R.ODU.1+rNC.2*R.ODU.2+rNC.3*R.ODU.3) -
        por.grid$mid/phi.N$mid*(R.NOx+R.ana) + alpha$mid*(NH4.w-NH4)

      dNO3 <- trans.NO3$dC - svf.grid$mid/por.grid$mid*0.8*(R.NO3.1+R.NO3.2+R.NO3.3) + R.NOx - R.ana + alpha$mid*(NO3.w-NO3)

      dODU <- trans.ODU$dC + svf.grid$mid/por.grid$mid*(R.ODU.1+R.ODU.2+R.ODU.3) - R.ODUox - svf.grid$mid/por.grid$mid*OduDepo + alpha$mid*(ODU.w-ODU)

      dDIC  <- trans.DIC$dC + svf.grid$mid/por.grid$mid*(R.DIC.1+R.DIC.2+R.DIC.3) + alpha$mid*(DIC.w-DIC)
      dTAlk <- trans.TAlk$dC + R.TAlk + alpha$mid*(TAlk.w-TAlk)

      # diagnostics: the v2 set (same formulas) + a few from v1's transient
      # fluxes in umol cm-2 yr-1 (x 10/365 -> mmol m-2 d-1)
      return(list(c(dOM1=dOM1, dOM2=dOM2, dOM3=dOM3,
                    dO2=dO2, dNH4=dNH4, dNO3=dNO3,
                    dODU=dODU, dDIC=dDIC, dTAlk=dTAlk),
                  o2_flux   = trans.O2$flux.up  + sum(por.grid$mid*grid$dx*alpha$mid*(O2.w-O2)),
                  nh4_flux  = trans.NH4$flux.up + sum(por.grid$mid*grid$dx*alpha$mid*(NH4.w-NH4)),
                  no3_flux  = trans.NO3$flux.up + sum(por.grid$mid*grid$dx*alpha$mid*(NO3.w-NO3)),
                  tot_denit = 0.8*sum(svf.grid$mid*grid$dx*(R.NO3.1+R.NO3.2+R.NO3.3)) + 2*sum(por.grid$mid*grid$dx*R.ana),
                  anammox   = 2*sum(por.grid$mid*grid$dx*R.ana),
                  dic_flux  = trans.DIC$flux.up + sum(por.grid$mid*grid$dx*alpha$mid*(DIC.w-DIC)),
                  OM_flux   = trans.OM1$flux.up+trans.OM2$flux.up+trans.OM3$flux.up,
                  OM_burial = trans.OM1$flux.down+trans.OM2$flux.down+trans.OM3$flux.down,
                  odu_flux  = trans.ODU$flux.up + sum(por.grid$mid*grid$dx*alpha$mid*(ODU.w-ODU)),
                  talk_flux = trans.TAlk$flux.up + sum(por.grid$mid*grid$dx*alpha$mid*(TAlk.w-TAlk)),
                  tot_remin = sum(svf.grid$mid*grid$dx*(R.DIC.1+R.DIC.2+R.DIC.3)),
                  nitrification = sum(por.grid$mid*grid$dx*R.NOx),
                  oxic_remin    = sum(svf.grid$mid*grid$dx*(R.O2.1+R.O2.2+R.O2.3)),
                  OC_storage    = sum((OM1+OM2+OM3)*svf.grid$mid*grid$dx)/1e6*1e4,
                  Db_top    = Db.grid$mid[1],
                  alpha_top = alpha$mid[1]
      ))
    })
  }

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
