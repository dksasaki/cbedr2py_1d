"""
Newton-Krylov spin-up and 100-yr drift check for the five field sites.

Same workflow as test_krylov.py; params and f_year are built per site from
the boundary-condition table (BW = bottom water). Unit conversions:
  concentrations  uM -> mM               (/1000)
  OC flux         mmol C m-2 d-1 -> umol C cm-2 yr-1   (x36.5)
Not in the table (taken from test_krylov.py): porosity, OM fractions, TAlk.w.
'depth' is not used by the transient model, it is kept only as metadata.
"""
import time
import numpy as np
import matplotlib.pyplot as plt
import cbed2py
from cbed2py.spinup import newton_krylov_spinup, periodic_forcing, annual_budget

cbed2py.set_model('fast')

FLUX_CONV = 36.5          # mmol C m-2 d-1 -> umol C cm-2 yr-1
MODE = 'instantaneous'    # J.OM mode for reactivity, bioturbation, bioirrigation


def site_params(J_mmol, temp, O2_uM, NO3_uM, NH4_uM, DIC_uM, w, depth=None):
    """CBED params (CBED units) for one site, from the table units."""
    return {
        'por.0': 0.75, 'por.inf': 0.75, 'depth': depth,
        'f.OM1': 0.33, 'f.OM2': 0.62, 'f.OM3': 0.05,
        'J.OM': J_mmol * FLUX_CONV,
        'O2.w': O2_uM / 1e3, 'NO3.w': NO3_uM / 1e3, 'NH4.w': NH4_uM / 1e3,
        'DIC.w': DIC_uM / 1e3, 'TAlk.w': 2.25,
        'temp': temp, 'S': 34.0, 'w': w,
        'jpoc_mode_reactivity':    MODE,
        'jpoc_mode_bioturbation':  MODE,
        'jpoc_mode_bioirrigation': MODE,
    }


#                       J.OM    T      O2    NO3   NH4   DIC   w
sites = {
    'site1': site_params(11.05,  6.87, 276,  8.0,  0.6,  2164, 0.03599, depth=74.0),
    'site2': site_params(2.78,  16.3,  137, 15.5,  0.5,  2250, 0.008067),
    'site3': site_params(10.8,  11.6,    0, 28.0,  0.8,  2358, 0.0337),
    'site4': site_params(0.12,   0.21, 214, 31.0,  0.05, 2344, 0.001390),
    'site5': site_params(0.2128, 1.73, 251, 21.0,  0.1,  2258, 0.003091),
}

# seasonal amplitude of bottom-water temperature [degC] (keep cold sites > 0 degC)
temp_amp = {'site1': 2.0, 'site2': 2.0, 'site3': 2.0, 'site4': 0.2, 'site5': 0.5}


# ----------------------------------------------------------------------------
site = 'site1'
# ----------------------------------------------------------------------------
params = sites[site]

# one forcing year (repeated by the spin-up)
t_year = np.arange(0, 1, 1/365)
f_year = {'time': t_year,
          'J.OM': params['J.OM'] * (1 + 0.5*np.sin(2*np.pi*t_year)),
          'temp': params['temp'] + temp_amp[site]*np.sin(2*np.pi*(t_year - 0.6))}

n_periods = 10 if site in ('site4', 'site5') else 1   # slow sites: longer shots


# spin-up to the periodic state (Newton-Krylov)
t0 = time.time()
y, hist, ds_last = newton_krylov_spinup(params, f_year, n_periods=n_periods, max_newton=10,
                                        checkpoint=f'spinup_nk_{site}_{MODE}.npz')
print(f'{site}: spin-up {hist[-1]["runs"] * n_periods} model years, {time.time() - t0:.0f} s')
if not hist[-1]['converged']:
    print(f'WARNING: not converged; restart with yini=np.load("spinup_nk_{site}_{MODE}.npz")["y"]')


# verification run: same seasonal forcing, should show no trend
t_exp = np.arange(0, 100 + 1e-9, 1/365)
f_exp = periodic_forcing(f_year, 0, 100)
ds = cbed2py.run_transient(t_exp, forcing=f_exp, params=params, yini=y)
# ds.to_netcdf(f'verify_{site}_{MODE}.nc')

imb = [annual_budget(ds, T)['imbalance'] for T in (1, 10, 50, 100)]
print('annual imbalance at yr 1, 10, 50, 100 [umol C cm-2 yr-1]:', np.round(imb, 4),
      f'({100*abs(imb[-1])/params["J.OM"]:.3f}% of J.OM)')
print('OC_storage change over 100 yr (relative):',
      float((ds.OC_storage[-1] - ds.OC_storage[0]) / ds.OC_storage[0]))


# % change of OM1-3 in the top, mid and bottom layers, sampled once a year
# at the same phase (t = 0, 1, 2, ...) so the seasonal cycle drops out
yrs = np.arange(0, 101)
ann = ds.sel(time=yrs, method='nearest')
N = ds.z.size
layers = {'top': 0, 'mid': N // 2, 'bottom': N - 1}

fig, axs = plt.subplots(1, 3, figsize=(11, 3.5), sharex=True, constrained_layout=True)
for ax, (name, k) in zip(axs, layers.items()):
    for s in ('OM1', 'OM2', 'OM3'):
        c = ann[s].isel(z=k)
        if c[0] > 1e-6:                      # skip pools that are ~0 there (e.g. OM1 at depth)
            ax.plot(yrs, 100 * (c - c[0]) / c[0], label=s)
    ax.axhline(0, color='0.6', lw=0.8)
    ax.set(title=f'{name} layer (z = {float(ds.z[k]):.2f} cm)', xlabel='time [yr]')
    ax.legend(frameon=False)
axs[0].set_ylabel('change vs year 0 [%]')
fig.suptitle(f'{site}: drift after spin-up')
# fig.savefig(f'drift_{site}_{MODE}.png', dpi=150)
plt.show()


# Hovmoller (time x depth) of the % change of OM1-3 vs year 0, same annual
# sampling. Each panel has its own symmetric colour range (drift sizes differ
# a lot between pools); layers where the pool is ~0 are masked.
fig, axs = plt.subplots(1, 3, figsize=(12, 4), sharey=True, constrained_layout=True)
for ax, s in zip(axs, ('OM1', 'OM2', 'OM3')):
    c = ann[s]                                    # (time, z)
    c0 = c.isel(time=0)
    pct = 100 * (c - c0) / c0.where(c0 > 1e-3 * c0.max())
    vmax = float(np.nanmax(np.abs(pct))) or 1e-12
    pm = ax.pcolormesh(yrs, ds.z, pct.T, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                       shading='nearest')
    fig.colorbar(pm, ax=ax, label='change vs year 0 [%]')
    ax.set(title=s, xlabel='time [yr]')
axs[0].invert_yaxis()                             # shared y: inverts all panels
axs[0].set_ylabel('depth [cm]')
fig.suptitle(f'{site}: OM drift after spin-up (same phase each year)')
# fig.savefig(f'hovmoller_{site}_{MODE}.png', dpi=150)
plt.show()


# Hovmoller (time x depth) of the OM1-3 concentrations, full output
# resolution (seasonal cycle included). Log colour scale: concentrations
# drop by orders of magnitude with depth; values ~0 are masked.
from matplotlib.colors import LogNorm

fig, axs = plt.subplots(3,1, figsize=(12, 4), sharey=True, constrained_layout=True)
for ax, s in zip(axs, ('OM1', 'OM2', 'OM3')):
    c = ds[s]                                     # (time, z)
    cmax = float(c.max())
    cm = c.where(c > 1e-6 * cmax)
    pm = ax.pcolormesh(ds.time, ds.z, cm.T, cmap='viridis', shading='nearest',
                       norm=LogNorm(vmin=float(cm.min()), vmax=cmax))
    fig.colorbar(pm, ax=ax, label=f'{s} [umol cm-3]')
    ax.set(title=s, xlabel='time [yr]')
axs[0].invert_yaxis()                             # shared y: inverts all panels
axs[0].set_ylabel('depth [cm]')
fig.suptitle(f'{site}: OM concentrations')
# fig.savefig(f'hovmoller_conc_{site}_{MODE}.png', dpi=150)
plt.show()



# benthic fluxes, converted umol cm-2 yr-1 -> mmol m-2 d-1 (/36.5)
# sign as in the model: positive into the sediment
fig, axs = plt.subplots(4, 1, figsize=(12, 6), sharex=True, constrained_layout=True)
for ax, s in zip(axs, ('o2_flux', 'dic_flux', 'nh4_flux', 'no3_flux')):
    ax.plot(ds.time, ds[s] / FLUX_CONV)
    ax.axhline(0, color='0.6', lw=0.8)
    ax.set(title=s, ylabel='mmol m$^{-2}$ d$^{-1}$')
axs[-1].set_xlabel('time [yr]')
fig.suptitle(f'{site}: benthic fluxes (positive into the sediment)')
plt.show()

