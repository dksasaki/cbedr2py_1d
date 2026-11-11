import numpy as np
import cbed2py
import glob
import xarray as xr

def save_experiment_config(m):
    fout = f'example_{m}.nc'

    if not glob.glob(fout):
        # point values and constant (spin-up) forcing, CBED units
        params = {
            'por.0': 0.75, 'por.inf': 0.75,          # porosity
            'depth': 250.0,                          # m
            'f.OM1': 0.33, 'f.OM2': 0.62, 'f.OM3': 0.05,
            'J.OM': 400.0, 'O2.w': 0.25, 'NO3.w': 0.015, 'NH4.w': 0.001,
            'DIC.w': 2.15, 'TAlk.w': 2.25, 'temp': 7.0, 'S': 33.0, 'w': 0.05,
            'jpoc_mode_reactivity':    m,   
            'jpoc_mode_bioturbation':  m,
            'jpoc_mode_bioirrigation': m,
            # 'J.avg.days': 7,                            # window for 'trailing_mean'
        }

        # time-varying forcing, CBED units, time in years
        t = np.arange(0, 100 + 1e-9, 1/365)
        forcing = {'time': t,
                'J.OM': 400 * (1 + 0.5*np.sin(2*np.pi*t)),
                'temp': 7 + 2*np.sin(2*np.pi*(t - 0.6))}

        ds = cbed2py.run_transient(t, forcing=forcing, params=params)   # spin-up from steady state
        ds.to_netcdf(fout)
    else:
        ds = xr.open_dataset(fout)
    return ds
if __name__ == '__main__':

    cbed2py.set_model('fast')

    mode_pars = ['instantaneous',
                'longterm_mean',
                'trailing_mean']


    exps = {}
    for m in mode_pars:
        exps[m] = save_experiment_config(m)

