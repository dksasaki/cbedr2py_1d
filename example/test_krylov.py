import time
import numpy as np, cbed2py
from cbed2py.spinup import newton_krylov_spinup, periodic_forcing, annual_budget

cbed2py.set_model('fast')


# one forcing year (repeated by the spin-up)
t_year = np.arange(0, 1, 1/365)
f_year = {'time': t_year,
          'J.OM': 400*(1 + 0.5*np.sin(2*np.pi*t_year)),
          'temp': 7 + 2*np.sin(2*np.pi*(t_year - 0.6))}


params = {
       'por.0': 0.75, 'por.inf': 0.75, 'depth': 250.0,
       'f.OM1': 0.33, 'f.OM2': 0.62, 'f.OM3': 0.05,
       'J.OM': 400.0, 'O2.w': 0.25, 'NO3.w': 0.015, 'NH4.w': 0.001,
       'DIC.w': 2.15, 'TAlk.w': 2.25, 'temp': 7.0, 'S': 33.0, 'w': 0.05,
       'jpoc_mode_reactivity':    'trailing_mean',
       'jpoc_mode_bioturbation':  'trailing_mean',
       'jpoc_mode_bioirrigation': 'trailing_mean',
   }


# spin-up to the periodic state (Newton-Krylov)
t0 = time.time()
y, hist, ds_last = newton_krylov_spinup(params, f_year, max_newton=8,
                                        checkpoint='spinup_nk_trailing_mean_fast.npz')
print(f'spin-up: {hist[-1]["runs"]} model years, {time.time() - t0:.0f} s')
if not hist[-1]['converged']:
    print('WARNING: spin-up not converged; restart with yini=np.load(...)["y"]')


# verification run: same seasonal forcing, should show no trend
t_exp = np.arange(0, 100 + 1e-9, 1/365)
f_exp = periodic_forcing(f_year, 0, 100)          # with lead-in for the trailing mean
ds = cbed2py.run_transient(t_exp, forcing=f_exp, params=params, yini=y)

imb = [annual_budget(ds, T)['imbalance'] for T in (1, 10, 50, 100)]
print('annual imbalance at yr 1, 10, 50, 100 [umol C cm-2 yr-1]:', np.round(imb, 4))
print('OC_storage change over 100 yr (relative):',
      float((ds.OC_storage.values[-1] - ds.OC_storage.values[0]) / ds.OC_storage.values[0]))
