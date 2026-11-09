"""
test_wrapper.py -- checks of the transient wrapper.

Run  pixi run test        (or: pixi run python tests/test_wrapper.py)

1. steady state of the transient equations == cbed_model() of cbed_v2_func.R
2. constant forcing, started from the steady state, does not drift
3. a sinusoidal J.OM forcing runs and the fluxes respond
4. changing a coefficient (k2.0) changes the solution
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
