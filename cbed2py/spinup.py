"""
cbed2py.spinup -- spin-up to the periodic (cyclostationary) state with
Aitken (delta-squared) extrapolation.

Theory
  With forcing of period P, the long-term state is a fixed point of the
  one-year map Phi: y0 -> y(P). Near it, the error after n years is
  e_n = sum_i c_i mu_i^n v_i (mu_i: Floquet multipliers). Fast porewater
  modes have mu ~ 0 and vanish within 1-2 years; the OM inventory has
  mu = exp(-P/tau) with tau ~ centuries, which makes brute-force spin-up slow.
  Once only the slow mode is left, any state variable sampled every dt years
  (at the same phase of the forcing) follows
        x_k = x* + A q^k,   q = mu^dt,
  and three samples x0, x1, x2 give
        q  = (x2 - x1) / (x1 - x0)
        x* = x2 + (x2 - x1) * q / (1 - q)      (= Aitken's delta-squared)
  Restarting from x* and repeating removes the remaining (other-mode) error.

Functions
  periodic_forcing(forcing, t_start, t_end, period=1, lead=None)
      repeat one period of forcing over [t_start - lead, t_end]
  aitken_extrapolate(y0, y1, y2, ...)
      the extrapolation itself (state vectors of 9*N)
  aitken_from_dataset(ds, t0, dt, ...)
      extrapolate from a run already done (e.g. a saved netCDF)
  aitken_spinup(params, forcing, yini=None, dt=20, ...)
      run -> extrapolate -> restart, until the drift is below tolerance
  annual_budget(ds, t_end, period=1)
      annual means of the OM budget terms over [t_end - period, t_end]
  newton_krylov_spinup(params, forcing, yini=None, ...)
      alternative: solve Phi(y) - y = 0 directly with Jacobian-free
      Newton-Krylov (Newton-GMRES), see its docstring

Only the OM pools (OM1-3) are extrapolated by default: the porewater
species are slaved to them and re-equilibrate within `skip` years of each
restart.
"""
import numpy as np
import xarray as xr

from .wrapper import run_transient, steady_state, STATES, FORCING

SLOW_STATES = ('OM1', 'OM2', 'OM3')


# ---------------------------------------------------------------------------
# forcing
# ---------------------------------------------------------------------------
def periodic_forcing(forcing, t_start, t_end, period=1.0, lead=None):
    """
    Repeat one period of forcing over [t_start - lead, t_end].

    forcing : dict / DataFrame with 'time' [yr] covering one period (time
              taken modulo `period`) and any FORCING columns.
    lead    : forcing kept before t_start, so that the trailing mean of
              J.OM (mode 'trailing_mean') sees periodic forcing from the
              first day. Default: one period.
    The sample spacing of the input is kept.
    """
    if lead is None:
        lead = period
    tp = np.asarray(forcing['time'], dtype=float) % period
    order = np.argsort(tp)
    tp, keep = np.unique(tp[order], return_index=True)
    dtf = np.median(np.diff(tp)) if tp.size > 1 else period / 365
    t = np.arange(t_start - lead, t_end + dtf / 2, dtf)
    out = {'time': t}
    for v in FORCING:
        if v in forcing:
            fp = np.asarray(forcing[v], dtype=float)[order][keep]
            out[v] = np.interp(t, tp, fp, period=period)
    return out


# ---------------------------------------------------------------------------
# state vector <-> named pieces
# ---------------------------------------------------------------------------
def _state_vector(ds_t):
    """state vector (9*N,) from a Dataset at one time"""
    return np.concatenate([ds_t[s].values for s in STATES])


def _slices(N):
    return {s: slice(k * N, (k + 1) * N) for k, s in enumerate(STATES)}


# ---------------------------------------------------------------------------
# Aitken
# ---------------------------------------------------------------------------
def aitken_extrapolate(y0, y1, y2, N, variables=SLOW_STATES, q_mode='species',
                       q_max=0.99, rel_min=1e-9):
    """
    Aitken extrapolation of state vectors sampled dt years apart at the same
    forcing phase.

    variables : states to extrapolate; the others are taken from y2.
    q_mode    : 'species'   one q per species, least squares over depth,
                            q = sum(d1*d2)/sum(d1^2)  (robust; assumes one
                            slow mode per species)
                'pointwise' one q per depth and species (classical Aitken)
    q_max     : q is accepted only in (0, q_max): monotone, convergent drift.
                Larger q (very slow mode vs dt) would make a large, uncertain
                jump; take a longer dt instead.
    rel_min   : increments |x1-x0| below rel_min*|x2| are treated as converged.

    Returns (y_star, info) with info[s] = dict(q=..., tau_over_dt=..., n_used=...).
    """
    y0, y1, y2 = (np.asarray(y, dtype=float) for y in (y0, y1, y2))
    ys = y2.copy()
    info = {}
    for s, sl in _slices(N).items():
        if s not in variables:
            continue
        x0, x1, x2 = y0[sl], y1[sl], y2[sl]
        d1, d2 = x1 - x0, x2 - x1
        active = np.abs(d1) > rel_min * np.maximum(np.abs(x2), 1e-30)

        if q_mode == 'species':
            den = np.sum(d1[active] ** 2)
            q = np.sum(d1[active] * d2[active]) / den if den > 0 else np.nan
            q = np.full(x0.shape, q)
        elif q_mode == 'pointwise':
            with np.errstate(divide='ignore', invalid='ignore'):
                q = np.where(active, d2 / d1, np.nan)
        else:
            raise ValueError("q_mode must be 'species' or 'pointwise'")

        ok = active & np.isfinite(q) & (q > 0) & (q < q_max)
        xs = x2.copy()
        xs[ok] = x2[ok] + d2[ok] * q[ok] / (1 - q[ok])
        xs = np.where(x2 >= 0, np.maximum(xs, 0.0), xs)   # keep concentrations >= 0
        ys[sl] = xs

        qv = q[ok]
        qm = float(np.median(qv)) if qv.size else np.nan
        info[s] = {'q': qm,
                   'tau_over_dt': -1 / np.log(qm) if np.isfinite(qm) else np.nan,
                   'n_used': int(ok.sum()), 'n': int(x0.size),
                   'max_jump_rel': float(np.max(np.abs(xs - x2) / np.maximum(np.abs(x2), 1e-30)))}
    return ys, info


def _at(ds, t):
    """Dataset at the output time nearest to t (warns if off by > 1 day)"""
    i = int(np.argmin(np.abs(ds.time.values - t)))
    if abs(ds.time.values[i] - t) > 1 / 365:
        print(f'WARNING: nearest output to t={t} is t={ds.time.values[i]:.4f}')
    return ds.isel(time=i)


def aitken_from_dataset(ds, t0, dt, **kw):
    """
    Extrapolate from a run already done: states at t0, t0+dt, t0+2dt
    (same forcing phase, i.e. dt a whole number of periods).
    Extra keywords go to aitken_extrapolate. Returns (y_star, info).
    """
    N = ds.sizes['z']
    ys = [_state_vector(_at(ds, t0 + k * dt)) for k in range(3)]
    y_star, info = aitken_extrapolate(*ys, N=N, **kw)
    for s in info:
        info[s]['tau_yr'] = info[s]['tau_over_dt'] * dt
    return y_star, info


# ---------------------------------------------------------------------------
# diagnostics
# ---------------------------------------------------------------------------
def annual_budget(ds, t_end, period=1.0):
    """
    Annual means (trapezoid over the output times in [t_end-period, t_end])
    of the OM budget: OM_flux, tot_remin, OM_burial and
    imbalance = OM_flux - tot_remin - OM_burial  [umol C cm-2 yr-1],
    which is d(OM inventory)/dt and is 0 in the periodic state.
    """
    t = ds.time.values
    m = (t >= t_end - period - 1e-9) & (t <= t_end + 1e-9)
    tt = t[m]
    if tt.size < 3:
        raise ValueError('need sub-annual output over the last period')

    def mean(v):
        return float(np.trapz(ds[v].values[m], tt) / (tt[-1] - tt[0]))

    b = {v: mean(v) for v in ('OM_flux', 'tot_remin', 'OM_burial')}
    b['imbalance'] = b['OM_flux'] - b['tot_remin'] - b['OM_burial']
    return b


# ---------------------------------------------------------------------------
# spin-up driver
# ---------------------------------------------------------------------------
def aitken_spinup(params, forcing, yini=None, dt=20, skip=2, n_rounds=6,
                  tol_imbalance=0.05, tol_drift=1e-6, period=1.0,
                  n_sub=12, verbose=True, checkpoint=None, **kw):
    """
    Spin up CBED to its periodic state under periodic forcing.

    Each round: from y, run skip + 2*dt periods with output n_sub times per
    period; take the states at t1 = skip, t1+dt, t1+2dt (whole periods: same
    forcing phase); extrapolate the OM pools with Aitken; restart.

    params        : as for run_transient
    forcing       : one period of forcing (dict with 'time' + FORCING columns),
                    repeated with periodic_forcing()
    yini          : start state (None -> steady_state(params))
    dt            : sample spacing [periods]. Choose dt ~ tau/5 .. tau/2 of the
                    slow mode (q = exp(-dt/tau) ~ 0.6 - 0.8); tau ~ 200 yr here
    skip          : periods run before the first sample, so the porewater
                    has re-equilibrated after the restart
    tol_imbalance : stop when |annual-mean OM_flux - tot_remin - OM_burial|
                    over the last period < this [umol C cm-2 yr-1]
    tol_drift     : ... or when max relative OM change per period < this
    checkpoint    : path of an .npz file; after every round it holds the
                    current state `y` (restart with yini=np.load(path)['y'])
                    and the per-round imbalance and drift
    kw            : passed to aitken_extrapolate (variables, q_mode, q_max)

    Returns (y, history, ds_last):
      y        state vector (9*N), usable as yini of run_transient
      history  list of dicts per round (budget, drift, q, tau, ...)
      ds_last  the last run (before the last extrapolation)
    """
    if yini is None:
        yini = steady_state(params)
    if isinstance(yini, xr.Dataset):
        yini = _state_vector(yini.isel(time=-1) if 'time' in yini.dims else yini)
    y = np.asarray(yini, dtype=float).ravel()

    t1 = float(skip) * period
    t_end = t1 + 2 * dt * period
    t_anchor = [t1, t1 + dt * period, t_end]
    # Output every period/n_sub over the whole run, not only at the anchors:
    # lsoda (ode.1D) allows at most maxsteps = 5000 internal steps between two
    # output times, and with daily-varying forcing a gap of decades exceeds it.
    n_out = int(round(t_end / period * n_sub))
    times = np.round(np.linspace(0.0, t_end, n_out + 1), 10)
    frc = periodic_forcing(forcing, 0.0, t_end, period=period)

    history = []
    ds = None
    for k in range(n_rounds):
        ds = run_transient(times, forcing=frc, params=params, yini=y)
        if ds.time.values[-1] < t_end - 1e-6:
            raise RuntimeError(
                f'round {k}: integration stopped at t={ds.time.values[-1]:.4g} '
                f'of {t_end:g} yr (see the R/deSolve messages above); '
                f'try a larger n_sub, or check the state y passed in')
        N = ds.sizes['z']
        ys = [_state_vector(_at(ds, t)) for t in t_anchor]
        sl = _slices(N)

        bud = annual_budget(ds, t_end, period)
        om = np.concatenate([ys[2][sl[s]] for s in SLOW_STATES])
        dom = np.concatenate([ys[2][sl[s]] - ys[1][sl[s]] for s in SLOW_STATES])
        drift = float(np.max(np.abs(dom) / np.maximum(np.abs(om), 1e-30)) / dt)

        rec = {'round': k, 'years_run': (k + 1) * t_end,
               'imbalance': bud['imbalance'], 'tot_remin': bud['tot_remin'],
               'OM_burial': bud['OM_burial'], 'drift_per_period': drift}

        if abs(bud['imbalance']) < tol_imbalance or drift < tol_drift:
            rec['converged'] = True
            history.append(rec)
            y = ys[2]
            _save_checkpoint(checkpoint, y, history)
            if verbose:
                _report(rec)
                print('converged')
            break

        y, info = aitken_extrapolate(*ys, N=N, **kw)
        rec['converged'] = False
        rec['aitken'] = {s: {**v, 'tau_yr': v['tau_over_dt'] * dt * period}
                         for s, v in info.items()}
        history.append(rec)
        _save_checkpoint(checkpoint, y, history)
        if verbose:
            _report(rec)
    else:
        if verbose:
            print(f'not converged after {n_rounds} rounds')
    return y, history, ds


def _save_checkpoint(path, y, history):
    """state + per-round budget, written after every round (if path is set)"""
    if path is None:
        return
    np.savez(path, y=y,
             imbalance=np.array([h['imbalance'] for h in history]),
             drift=np.array([h['drift_per_period'] for h in history]),
             converged=history[-1]['converged'])


def _report(rec):
    s = (f"round {rec['round']}: imbalance {rec['imbalance']:+.4f} umol cm-2 yr-1, "
         f"drift {rec['drift_per_period']:.2e} /period")
    if 'aitken' in rec:
        s += '  ' + ', '.join(f"{k}: q={v['q']:.3f} tau={v['tau_yr']:.0f}yr"
                              f" jump={v['max_jump_rel']:.1e}"
                              for k, v in rec['aitken'].items())
    print(s, flush=True)


# ---------------------------------------------------------------------------
# Newton-Krylov spin-up
# ---------------------------------------------------------------------------
def _gmres(Av, b, m, rtol):
    """
    GMRES (no restart, modified Gram-Schmidt, Givens rotations) for A x = b
    with A given only as products Av(v). Stops after m iterations or when
    ||b - A x|| <= rtol*||b||. Returns (x, relative residual, iterations).
    """
    n = b.size
    beta = np.linalg.norm(b)
    if beta == 0:
        return np.zeros(n), 0.0, 0
    V = np.zeros((m + 1, n))
    H = np.zeros((m + 1, m))
    cs, sn = np.zeros(m), np.zeros(m)
    g = np.zeros(m + 1)
    g[0] = beta
    V[0] = b / beta
    k = 0
    for j in range(m):
        w = Av(V[j])
        for i in range(j + 1):
            H[i, j] = w @ V[i]
            w = w - H[i, j] * V[i]
        H[j + 1, j] = np.linalg.norm(w)
        if H[j + 1, j] > 0:
            V[j + 1] = w / H[j + 1, j]
        for i in range(j):                                 # previous rotations
            h = cs[i] * H[i, j] + sn[i] * H[i + 1, j]
            H[i + 1, j] = -sn[i] * H[i, j] + cs[i] * H[i + 1, j]
            H[i, j] = h
        r = np.hypot(H[j, j], H[j + 1, j])                 # new rotation
        cs[j], sn[j] = (1.0, 0.0) if r == 0 else (H[j, j] / r, H[j + 1, j] / r)
        H[j, j], H[j + 1, j] = r, 0.0
        g[j + 1] = -sn[j] * g[j]
        g[j] = cs[j] * g[j]
        k = j + 1
        if abs(g[k]) <= rtol * beta or H[j, j] == 0:
            break
    yk = np.linalg.solve(np.triu(H[:k, :k]), g[:k]) if k else np.zeros(0)
    return V[:k].T @ yk, abs(g[k]) / beta, k


def newton_krylov_spinup(params, forcing, yini=None, period=1.0, n_periods=1,
                         n_sub=12, max_newton=10, gmres_iter=20, eta=0.1,
                         fd_rel=1e-3, scale_floor=1e-3, tol_imbalance=0.05,
                         tol_f=3e-5, max_backtrack=4, checkpoint=None,
                         verbose=True):
    """
    Spin up CBED to its periodic state with Jacobian-free Newton-Krylov.

    Theory
      The periodic state is the root of G(y) = Phi(y) - y, Phi = the map over
      n_periods forcing periods. Newton: solve J_G d = -G, J_G = M - I
      (M: monodromy matrix), with GMRES. J_G is never formed: each GMRES
      iteration needs one product, by finite differences,
            J_G v ~ (G(y + e v) - G(y)) / e      (one model run of n_periods).
      Fast porewater modes give eigenvalues of J_G near -1, the slow OM modes
      eigenvalues near 0 (= mu - 1); GMRES needs ~ one iteration per
      eigenvalue cluster, so a few to ~20 runs per Newton step, and Newton
      removes the slow modes in one step instead of following them for
      centuries.

    Implementation
      - scaled unknowns x = y / s, s = max(|y_ini|, scale_floor * max|y_ini|
        of that species): every unknown is O(1), so GMRES and the norms are
        not dominated by OM2 (~1e3-1e4) or blind to O2 (~0 at depth).
      - finite-difference step: RMS perturbation of x = fd_rel. It must be
        well above the solver noise of Phi (lsoda rtol = atol = 1e-6), so
        fd_rel ~ 1e-4 - 1e-3; too small gives noisy J v, too large adds
        nonlinearity error (both only slow Newton, they do not bias the root).
      - inexact Newton: GMRES stops at relative residual eta.
      - backtracking line search on ||G||; concentrations kept >= 0.

    params, forcing : as for aitken_spinup (forcing = one period)
    yini          : start state (None -> steady_state(params)); a state that
                    is already close (e.g. from aitken_spinup, or a few
                    decades of run) needs fewer Newton steps
    n_periods     : periods per shot (1 is usual)
    n_sub         : output times per period (also limits lsoda's steps
                    between outputs, see aitken_spinup)
    max_newton    : maximum Newton steps
    gmres_iter    : maximum GMRES iterations (model runs) per Newton step
    tol_imbalance : stop when |annual-mean OM_flux - tot_remin - OM_burial|
                    over the last shot < this [umol C cm-2 yr-1] ...
    tol_f         : ... and RMS of the scaled residual G(y)/s per period < this.
                    Keep it above the solver-noise floor of G: with lsoda's
                    default rtol = atol = 1e-6, |G|/|y| ~ 1e-6 - 1e-5 even at
                    the periodic state, so tol_f < ~1e-5 is never reached
    checkpoint    : .npz path, written after every Newton step (y, history)

    Returns (y, history, ds_last): state vector usable as yini, per-step
    records, and the last one-shot run (from y).
    """
    if yini is None:
        yini = steady_state(params)
    if isinstance(yini, xr.Dataset):
        yini = _state_vector(yini.isel(time=-1) if 'time' in yini.dims else yini)
    y0 = np.asarray(yini, dtype=float).ravel()

    t_end = n_periods * period
    times = np.round(np.linspace(0.0, t_end, int(round(n_periods * n_sub)) + 1), 10)
    frc = periodic_forcing(forcing, 0.0, t_end, period=period)

    n = y0.size
    N = n // len(STATES)
    s = np.empty(n)
    for sl in _slices(N).values():
        ys = np.abs(y0[sl])
        s[sl] = np.maximum(ys, scale_floor * max(ys.max(), 1e-30))
    nonneg = y0 >= 0
    n_runs = [0]

    def shoot(x):
        """scaled residual G and the run, from scaled state x"""
        y = s * x
        ds = run_transient(times, forcing=frc, params=params, yini=y)
        n_runs[0] += 1
        if ds.time.values[-1] < t_end - 1e-6:
            raise RuntimeError(f'integration stopped at t={ds.time.values[-1]:.4g} '
                               f'of {t_end:g} yr (see the R/deSolve messages above)')
        return (_state_vector(ds.isel(time=-1)) - y) / s, ds

    def clip(x):
        return np.where(nonneg, np.maximum(x, 0.0), x)

    def rms(v):
        return float(np.sqrt(np.mean(v ** 2)))

    x = y0 / s
    F, ds = shoot(x)
    history = []
    for it in range(max_newton + 1):
        bud = annual_budget(ds, t_end, period)
        rec = {'newton': it, 'runs': n_runs[0], 'imbalance': bud['imbalance'],
               'rms_G': rms(F) / n_periods, 'max_G': float(np.max(np.abs(F))) / n_periods}
        converged = abs(bud['imbalance']) < tol_imbalance and rec['rms_G'] < tol_f
        rec['converged'] = converged
        if verbose:
            print(f"newton {it}: imbalance {rec['imbalance']:+.4f} umol cm-2 yr-1, "
                  f"|G|rms {rec['rms_G']:.2e}, |G|max {rec['max_G']:.2e} /period, "
                  f"model runs {n_runs[0]}", flush=True)
        if converged or it == max_newton:
            history.append(rec)
            _save_checkpoint(checkpoint, s * x, [{**h, 'drift_per_period': h['rms_G']} for h in history])
            if verbose:
                print('converged' if converged else f'not converged after {max_newton} Newton steps')
            break

        # Newton direction: J_G d = -G by GMRES with finite-difference products
        def Av(v):
            nv = np.linalg.norm(v)
            if nv == 0:
                return np.zeros_like(v)
            e = fd_rel * np.sqrt(n) / nv
            Fe, _ = shoot(x + e * v)
            return (Fe - F) / e

        d, relres, k = _gmres(Av, -F, gmres_iter, eta)
        rec.update(gmres_iter=k, gmres_relres=relres)

        # backtracking line search on ||G||
        lam, best = 1.0, None
        for _ in range(max_backtrack + 1):
            xn = clip(x + lam * d)
            Fn, dsn = shoot(xn)
            if best is None or np.linalg.norm(Fn) < np.linalg.norm(best[1]):
                best = (xn, Fn, dsn, lam)
            if np.linalg.norm(Fn) <= (1 - 1e-4 * lam) * np.linalg.norm(F):
                break
            lam /= 2
        else:
            if verbose:
                print(f'  line search: no sufficient decrease, taking lambda={best[3]:g}')
        x, F, ds, lam = best
        rec['step'] = lam
        history.append(rec)
        _save_checkpoint(checkpoint, s * x, [{**h, 'drift_per_period': h['rms_G']} for h in history])
        if verbose:
            print(f'  GMRES {k} iterations (rel. residual {relres:.1e}), step {lam:g}', flush=True)

    return s * x, history, ds
