"""
test_wrapper.py -- checks of the transient wrapper.

Run  pixi run test        (or: pixi run python tests/test_wrapper.py)

1. steady state of the transient equations == cbed_model() of cbed_v2_func.R
2. constant forcing, started from the steady state, does not drift
3. a sinusoidal J.OM forcing runs and the fluxes respond
4. changing a coefficient (k2.0) changes the solution
5. J modes: the J seen by k / Db / irrigation is J.OM(t), the long-term
   mean, or the trailing mean over J.avg.days (checked against numpy)
"""
import os.path as osp
import time

import numpy as np

from cbed2py import wrapper as cw      # module, to reach the R helpers (_r, _to_py)

r = cw._r()

# 1. ------------------------------------------------------------------------
r(f'source("{osp.join(cw.MODEL_DIR, "cbed_v2_func.R")}")')
r('v2 <- cbed_model(get_default_pars())')
y_v2 = cw._to_py(r('v2$y'))                                   # (N, 9)
ss   = cw.steady_state()
print('1. steady state, transient equations vs cbed_v2_func.R')
for k, s in enumerate(cw.STATES):
    ref = y_v2[:, k]
    rel = np.max(np.abs(ss[s].values - ref)) / max(np.max(np.abs(ref)), 1e-30)
    print(f'   {s:5s} max rel diff = {rel:.2e}')
for f in ('o2_flux', 'dic_flux', 'OM_burial'):
    print(f'   {f:10s} v2 = {float(cw._to_py(r(f"v2${f}"))[0]):.6g}   new = {float(ss[f]):.6g}')

# 2. ------------------------------------------------------------------------
t = np.arange(0, 1 + 1e-9, 7/365)
t0 = time.time()
out = cw.run_transient(t, forcing=None, yini=ss)
print(f'2. constant forcing, 1 yr ({time.time()-t0:.1f} s)')
for s in cw.STATES:
    ref = ss[s].values
    rel = np.max(np.abs(out[s].values - ref)) / max(np.max(np.abs(ref)), 1e-30)
    print(f'   {s:5s} max rel drift = {rel:.2e}')

# 3. ------------------------------------------------------------------------
p   = cw.default_params()
t   = np.arange(0, 2 + 1e-9, 1/365)
frc = {'time': t, 'J.OM': p['J.OM'] * (1 + 0.5*np.sin(2*np.pi*t))}
t0 = time.time()
out = cw.run_transient(t, forcing=frc, yini=ss)
print(f'3. sinusoidal J.OM, 2 yr daily ({time.time()-t0:.1f} s), {out.time.size} outputs')
for f in ('forcing_J.OM', 'OM_flux', 'o2_flux', 'dic_flux', 'tot_denit'):
    v = out[f].sel(time=slice(1, 2))
    print(f'   {f:13s} yr 2: min = {float(v.min()):9.4g}  max = {float(v.max()):9.4g}')

# 4. ------------------------------------------------------------------------
ss2 = cw.steady_state({'k2.0': 2 * p['k2.0']})
print('4. k2.0 x 2:  OM2 surface', f'{float(ss.OM2[0]):.4g} -> {float(ss2.OM2[0]):.4g}',
      f'| o2_flux {float(ss.o2_flux):.4g} -> {float(ss2.o2_flux):.4g}')

# 5. ------------------------------------------------------------------------
t5   = np.arange(0, 0.25 + 1e-9, 1/365)
frc5 = {'time': t5, 'J.OM': p['J.OM'] * (1 + 0.5*np.sin(2*np.pi*12*t5))}  # monthly cycle
par5 = {'jpoc_mode_reactivity': 'trailing_mean', 'jpoc_mode_bioturbation': 'longterm_mean',
        'jpoc_mode_bioirrigation': 'instantaneous', 'J.avg.days': 7}
t0 = time.time()
out5 = cw.run_transient(t5, forcing=frc5, params=par5, yini=ss)
# reference trailing mean: average of the interpolated forcing on a fine grid
tau  = 7/365
ref  = np.array([np.interp(np.linspace(s - tau, s, 2001), t5, frc5['J.OM']).mean()
                 for s in out5.time.values])
J_inst = np.interp(out5.time.values, t5, frc5['J.OM'])
print(f'5. J modes, 3 months daily ({time.time()-t0:.1f} s)')
print(f'   J_k   (trailing_mean) max rel err vs numpy = {np.max(np.abs(out5.J_k.values/ref - 1)):.1e}')
print(f'   J_Db  (longterm_mean) max |J_Db - J.OM|    = {np.max(np.abs(out5.J_Db.values - p["J.OM"])):.1e}')
print(f'   J_irr (instantaneous) max |J_irr - J.OM(t)| = {np.max(np.abs(out5.J_irr.values - J_inst)):.1e}')
print(f'   J.OM(t) range {J_inst.min():.0f}-{J_inst.max():.0f}, trailing mean range '
      f'{out5.J_k.min().item():.0f}-{out5.J_k.max().item():.0f}')
