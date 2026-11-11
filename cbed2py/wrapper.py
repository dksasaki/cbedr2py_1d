"""
cbed2py.wrapper -- Python wrapper for the transient CBED model (v2 coefficients).

The model is R code in cbed2py/model/cbed_v2_transient.R; it is called
through rpy2. The functions below are also available as cbed2py.<name>.

Functions
  default_params()                      all parameters and coefficients with their defaults
  steady_state(params)                  steady state for constant forcing (initial condition)
  run_transient(times, forcing, params, yini)
                                        transient run -> xarray.Dataset

Example
  import numpy as np, cbed2py
  t   = np.arange(0, 3 + 1e-9, 1/365)                       # [yr], daily, 3 years
  frc = {'time': t, 'J.OM': 400*(1 + 0.5*np.sin(2*np.pi*t))}
  out = cbed2py.run_transient(t, forcing=frc,
                         params={'J.OM': 400, 'k2.0': 2e-3}) # spin-up with J.OM = 400

Units are those of the CBED code: depth z [cm], time [yr], solutes [mM],
J.OM [umol C cm-2 yr-1], w [cm yr-1], fluxes [umol cm-2 yr-1].
"""
import os.path as osp

import numpy as np
import pandas as pd
import xarray as xr
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter

HERE      = osp.dirname(osp.abspath(__file__))
MODEL_DIR = osp.join(HERE, 'model')

# Which POC flux J each J.OM-dependent process uses. Python name -> R name
# in model/cbed_v2_transient.R; the value is one of JPOC_MODES (a string
# here, an integer code in R).
JPOC_MODE_PARS = {
    'jpoc_mode_reactivity':    'J.mode.k',     # k1-k3
    'jpoc_mode_bioturbation':  'J.mode.Db',    # Db
    'jpoc_mode_bioirrigation': 'J.mode.irr',   # I0 (bioirrigation)
}
JPOC_MODES = {
    'instantaneous': 0,   # the forcing J.OM(t)
    'longterm_mean': 1,   # the constant params['J.OM']
    'trailing_mean': 2,   # mean of J.OM(t) over the previous params['J.avg.days'] days
}

STATES  = ["OM1", "OM2", "OM3", "O2", "NH4", "NO3", "ODU", "DIC", "TAlk"]
FORCING = ["J.OM", "O2.w", "NO3.w", "NH4.w", "DIC.w", "TAlk.w", "temp", "S", "w"]

UNITS = {
    'OM1': 'umol cm-3 (solid)', 'OM2': 'umol cm-3 (solid)', 'OM3': 'umol cm-3 (solid)',
    'O2': 'mM', 'NH4': 'mM', 'NO3': 'mM', 'ODU': 'mM', 'DIC': 'mM', 'TAlk': 'mM',
    'J.OM': 'umol C cm-2 yr-1', 'O2.w': 'mM', 'NO3.w': 'mM', 'NH4.w': 'mM',
    'DIC.w': 'mM', 'TAlk.w': 'mM', 'temp': 'degC', 'S': 'psu', 'w': 'cm yr-1',
}

# ---------------------------------------------------------------------------
# R session (sourced once, on first use)
# ---------------------------------------------------------------------------
_r_ready = False

# model versions: name -> R file in MODEL_DIR (same function names in each)
MODELS = {
    'v2':      'cbed_v2_transient.R',          # reference implementation (ReacTran calls)
    'fast':    'cbed_v2_transient_fast.R',     # same equations, optimized (~19x faster)
}
_model = 'v2'


def set_model(name):
    """
    Select the R model used by steady_state() / run_transient():
      'v2'       cbed_v2_transient.R (default): reference implementation,
                 written with ReacTran (tran.1D, setup.prop.1D)
      'fast'     cbed_v2_transient_fast.R: the same equations, optimized
                 (precomputed splines, inline transport); ~19x faster
    Both use the generic_CBED.F90 coefficients (diffusion coefficients,
    anammox stoichiometry); their right-hand sides agree to rounding.
    """
    global _model, _r_ready
    if name not in MODELS:
        raise ValueError(f'model must be one of {list(MODELS)}, got {name!r}')
    if name != _model:
        _model = name
        _r_ready = False            # re-source on next use
    return _model


def get_model():
    """name of the R model in use (see set_model)"""
    return _model


def _r():
    global _r_ready
    if not _r_ready:
        ro.r(f'source("{osp.join(MODEL_DIR, "default_CBED_pars.R")}")')
        ro.r(f'source("{osp.join(MODEL_DIR, MODELS[_model])}")')
        _r_ready = True
    return ro.r


def _to_py(robj):
    """R object -> numpy array (matrices keep their shape)."""
    with localconverter(ro.default_converter + numpy2ri.converter) as cv:
        return np.asarray(cv.rpy2py(robj))


def _r_list(p):
    """Python parameters -> R named list.
    J modes are renamed to their R names and passed as integer codes."""
    out = {}
    for k, v in p.items():
        if k in JPOC_MODE_PARS:
            out[JPOC_MODE_PARS[k]] = ro.IntVector([JPOC_MODES[v]])
        else:
            out[k] = ro.FloatVector(np.atleast_1d(np.asarray(v, dtype=float)))
    return ro.ListVector(out)


def _r_dataframe(forcing):
    """dict / DataFrame with 'time' + forcing columns -> R data.frame."""
    df = pd.DataFrame(forcing)
    if 'time' not in df:
        raise KeyError("forcing needs a 'time' column [yr]")
    unknown = set(df.columns) - set(FORCING) - {'time'}
    if unknown:
        raise KeyError(f'unknown forcing columns {sorted(unknown)}; allowed: {FORCING}')
    if np.any(np.diff(df['time'].values) <= 0):
        raise ValueError("forcing 'time' must be strictly increasing")
    df = df.astype(float)
    return ro.DataFrame({c: ro.FloatVector(df[c].values) for c in df.columns})


def _r_to_dict(rlist):
    out = {}
    for name, val in zip(rlist.names, rlist):
        arr = np.asarray(val, dtype=float)
        out[name] = arr.item() if arr.size == 1 else arr
    return out

# ---------------------------------------------------------------------------
# parameters
# ---------------------------------------------------------------------------
def default_params():
    """
    All parameters with their default values, as a dict:
      - get_default_pars()      (model/default_CBED_pars.R): forcing values,
        porosity, grid, half-saturations, ...
      - cbed_t_default_coefs()  (model/cbed_v2_transient.R): the v2
        constants (k1.0, k2.0, k3.0, k.exp, k.adj.*, Db.*, irr.xL,
        odu.bur.*, J.avg.days)
      - the POC flux J used by each J.OM-dependent process:
          jpoc_mode_reactivity     k1-k3
          jpoc_mode_bioturbation   Db
          jpoc_mode_bioirrigation  I0
        each one of
          'instantaneous'  the forcing J.OM(t)               (default)
          'longterm_mean'  the constant params['J.OM']
          'trailing_mean'  mean of J.OM(t) over the previous
                           params['J.avg.days'] days (default 7)
    """
    r = _r()
    p = _r_to_dict(r('get_default_pars()'))
    p.update(_r_to_dict(r('cbed_t_default_coefs()')))
    code2mode = {c: m for m, c in JPOC_MODES.items()}
    for py_name, r_name in JPOC_MODE_PARS.items():   # R code -> Python string
        p[py_name] = code2mode[int(p.pop(r_name))]
    return p


def _full_params(params):
    p = default_params()
    if params:
        unknown = set(params) - set(p)
        if unknown:
            raise KeyError(f'unknown parameters {sorted(unknown)}')
        for f in JPOC_MODE_PARS:
            if f in params and params[f] not in JPOC_MODES:
                raise ValueError(f'{f} must be one of {list(JPOC_MODES)}, got {params[f]!r}')
        p.update(params)
    return p

# ---------------------------------------------------------------------------
# steady state
# ---------------------------------------------------------------------------
def steady_state(params=None):
    """
    Steady state of the transient equations with constant forcing taken from
    `params` (dict of overrides on top of default_params()).

    Returns an xarray.Dataset with the 9 state profiles on z [cm] and the
    flux diagnostics as scalars.
    """
    r = _r()
    p = _full_params(params)
    ro.globalenv['py_pars'] = _r_list(p)
    r('py_ss <- cbed_model_ss_t(py_pars)')

    y = _to_py(r('py_ss$y'))                       # (N, 9)
    z = _to_py(r('py_ss$z'))
    ds = xr.Dataset({s: ('z', y[:, k]) for k, s in enumerate(STATES)},
                    coords={'z': z})
    for name in r('names(py_ss)'):
        if name not in ('y', 'z'):
            ds[name] = float(_to_py(r(f'py_ss[["{name}"]]'))[0])
    ds.attrs['steady'] = int(bool(r('attr(py_ss, "steady")')[0]))
    if not ds.attrs['steady']:
        print('WARNING: steady state not reached (runsteady)')
    _add_attrs(ds, p)
    return ds

# ---------------------------------------------------------------------------
# transient run
# ---------------------------------------------------------------------------
def run_transient(times, forcing=None, params=None, yini=None):
    """
    Run CBED in transient mode.

    times   : output times [yr], increasing (the first is the initial time)
    forcing : dict or DataFrame with 'time' [yr] and any of
              J.OM, O2.w, NO3.w, NH4.w, DIC.w, TAlk.w, temp, S, w (CBED units).
              Linearly interpolated in time, held constant past the ends.
              Variables not given stay constant at their `params` value.
              None -> constant forcing.
    params  : dict of overrides on top of default_params(): point values
              (por.0, por.inf, depth, f.OM1-3, ...) and coefficients
              (k1.0, k2.0, ...). The forcing values here (J.OM, O2.w, ...)
              are the constant values used for the spin-up.
    yini    : initial state. None -> steady_state(params).
              Or a Dataset from steady_state() / run_transient() (the last
              time is used, so a run can be continued), or an array of 9*N.

    Returns an xarray.Dataset with
      states (time, z), diagnostics (time), the forcing at the output
      times (time), and the parameters as attributes.
    """
    r = _r()
    p = _full_params(params)
    times = np.asarray(times, dtype=float)

    # initial condition
    if yini is None:
        yini = steady_state(p)
    if isinstance(yini, xr.Dataset):
        if 'time' in yini.dims:
            yini = yini.isel(time=-1)
        yini = np.concatenate([yini[s].values for s in STATES])
    yini = np.asarray(yini, dtype=float).ravel()
    if yini.size != 9 * int(p['N']):
        raise ValueError(f'yini has {yini.size} values, expected 9*N = {9*int(p["N"])}')

    ro.globalenv['py_pars']  = _r_list(p)
    ro.globalenv['py_times'] = ro.FloatVector(times)
    ro.globalenv['py_yini']  = ro.FloatVector(yini)
    ro.globalenv['py_forc']  = ro.NULL if forcing is None else _r_dataframe(forcing)
    r('py_tr <- cbed_model_t(py_pars, py_times, py_yini, py_forc)')

    t = _to_py(r('py_tr$time'))
    z = _to_py(r('py_tr$z'))
    ds = xr.Dataset(coords={'time': t, 'z': z})
    for s in STATES:
        ds[s] = (('time', 'z'), _to_py(r(f'py_tr[["{s}"]]')))
    for name in r('names(py_tr$diag)'):
        ds[name] = ('time', _to_py(r(f'py_tr$diag[["{name}"]]')))

    # forcing actually seen by the model, at the output times
    fdf = None if forcing is None else pd.DataFrame(forcing)
    for v in FORCING:
        if fdf is not None and v in fdf:
            ds[f'forcing_{v}'] = ('time', np.interp(t, fdf['time'].values, fdf[v].values))
        else:
            ds[f'forcing_{v}'] = ('time', np.full(t.size, float(p[v])))

    if t.size < times.size:
        print(f'WARNING: integration stopped at t={t[-1]:.4g} yr '
              f'({t.size} of {times.size} output times)')
    _add_attrs(ds, p)
    return ds


def _add_attrs(ds, p):
    ds['z'].attrs = {'units': 'cm', 'long_name': 'depth below sediment surface'}
    if 'time' in ds.coords:
        ds['time'].attrs = {'units': 'yr', 'long_name': 'model time'}
    for v in ds.data_vars:
        base = v.replace('forcing_', '')
        if base in UNITS:
            ds[v].attrs['units'] = UNITS[base]
        elif v.endswith('_flux') or v in ('tot_denit', 'anammox', 'OM_burial', 'tot_remin',
                                         'nitrification', 'oxic_remin'):
            ds[v].attrs['units'] = 'umol cm-2 yr-1'
    ds.attrs['model'] = _model
    ds.attrs['model_file'] = MODELS[_model]
    for k, v in p.items():
        ds.attrs[f'par_{k}'] = v
