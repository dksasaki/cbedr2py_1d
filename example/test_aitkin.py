import numpy as np, cbed2py
from cbed2py.spinup import aitken_spinup

cbed2py.set_model('fast')


t1 = np.arange(0, 1, 1/365)                                  # one forcing year
f1 = {'time': t1,
      'J.OM': 400*(1 + 0.5*np.sin(2*np.pi*t1)),
      'temp': 7 + 2*np.sin(2*np.pi*(t1 - 0.6))}


params = {
       'por.0': 0.75, 'por.inf': 0.75, 'depth': 250.0,
       'f.OM1': 0.33, 'f.OM2': 0.62, 'f.OM3': 0.05,
       'J.OM': 400.0, 'O2.w': 0.25, 'NO3.w': 0.015, 'NH4.w': 0.001,
       'DIC.w': 2.15, 'TAlk.w': 2.25, 'temp': 7.0, 'S': 33.0, 'w': 0.05,
       'jpoc_mode_reactivity':    'trailing_mean',
       'jpoc_mode_bioturbation':  'trailing_mean',
       'jpoc_mode_bioirrigation': 'trailing_mean',
   }



y, hist, ds_last = aitken_spinup(params, f1, dt=30, skip=2, n_rounds=6)

t1 = np.arange(0, 100, 1/365)                                  # one forcing year
f1 = {'time': t1,
      'J.OM': 400*(1 + 0.5*np.sin(2*np.pi*t1)),
      'temp': 7 + 2*np.sin(2*np.pi*(t1 - 0.6))}

ds = cbed2py.run_transient(t1, forcing=f1, params=params, yini=y)  # your experiment