"""
Newton-Krylov spin-up from a saved point forcing, using monthly means.

Same inputs as spinup_from_forcing.py (files written by extract_forcing.ipynb):
    <name>_params.json   params dict
    <name>_forcing.npz   f_year dict (time + J.OM, O2.w, ...), daily

Why monthly: the R model interpolates the forcing linearly, so every forcing
sample is a kink in the boundary conditions and lsoda shortens its step at
each one. Noisy daily bottom-water forcing costs ~7.6 s per model year;
12 monthly means cost ~0.6 s. Daily noise in the bottom water has little
effect on the slow OM pools that the spin-up solves for, so the spin-up is
done with monthly forcing and then checked with the daily forcing.

Steps
  1. daily forcing -> 12 monthly means (times at mid-month)
  2. Newton-Krylov spin-up with the monthly forcing
  3. check run: --verify-years years with the ORIGINAL daily forcing,
     starting from the monthly periodic state; prints the annual OM
     imbalance and the year-to-year OC_storage change; figure of the
     forcing (daily vs monthly) and the fluxes
  4. drift run: --drift-years years (default 100) from the spin-up state,
     with the ORIGINAL daily forcing (~8 s/yr, ~13 min for 100 yr). Any
     trend combines the residual drift of the spin-up with the small shift
     between the monthly and the daily periodic states. --drift-forcing
     monthly re-runs the spin-up's own periodic problem instead (~0.6 s/yr,
     pure spin-up drift). Same diagnostics and plots as ../test_krylov_sites.py

Usage (from this folder):
    python spinup_monthly.py forcing_ix294_iy100
    python spinup_monthly.py forcing_ix294_iy100 --verify-years 0     # skip check
    python spinup_monthly.py forcing_ix294_iy100 --n-periods 10 --max-newton 6
    python spinup_monthly.py forcing_ix294_iy100 --tol-imb 0.005     # 0.5% of J.OM
    python spinup_monthly.py forcing_ix294_iy100 --drift-years 0     # skip drift run

Writes <name>_spinup_monthly.npz (y, in cbed2py state order, and the
monthly forcing used), a checkpoint <name>_spinup_monthly_ckpt.npz and,
figures <name>_transient.png (check run), <name>_drift_layers.png,
<name>_drift_hovmoller.png, <name>_om_hovmoller.png, <name>_drift_fluxes.png.
"""
import argparse
import json
import time
import numpy as np
import cbed2py
from cbed2py.spinup import newton_krylov_spinup, periodic_forcing, annual_budget

ap = argparse.ArgumentParser()
ap.add_argument('name', help='file stem, e.g. forcing_ix294_iy100')
ap.add_argument('--model', default='fast', help="'fast' or 'v2'")
ap.add_argument('--n-periods', type=int, default=5)
ap.add_argument('--max-newton', type=int, default=10)
ap.add_argument('--gmres-iter', type=int, default=20)
ap.add_argument('--tol-f', type=float, default=1e-4,
                help='RMS relative state change per year at convergence; '
                     'slow-mode drift over 100 yr is at most ~100*tol_f (1e-4 -> ~1%%)')
ap.add_argument('--tol-imb', type=float, default=0.02,
                help='OM imbalance tolerance as a fraction of the mean J.OM (0.02 = 2%%)')
ap.add_argument('--verify-years', type=int, default=3,
                help='years of daily-forced check run after the spin-up (0 = skip)')
ap.add_argument('--drift-years', type=int, default=100,
                help='years of the drift run after the spin-up (0 = skip)')
ap.add_argument('--drift-forcing', default='daily', choices=('daily', 'monthly'),
                help='forcing of the drift run: daily (the original forcing, default) '
                     'or monthly (the spin-up forcing, ~13x faster)')
args = ap.parse_args()


# ---------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------
with open(f'{args.name}_params.json') as f:
    params = json.load(f)
d = np.load(f'{args.name}_forcing.npz')
f_daily = {k: d[k] for k in d.files}

# TAlk.w must be in mM (~2.3). Values in the thousands are umol/kg left
# unconverted by the extraction; fix them here (better: fix the notebook).
if params.get('TAlk.w', 0) > 100:
    print(f'NOTE: TAlk.w = {params["TAlk.w"]:.1f} looks like umol/kg; dividing by 1000 -> mM')
    params['TAlk.w'] /= 1e3
    if 'TAlk.w' in f_daily:
        f_daily['TAlk.w'] = f_daily['TAlk.w'] / 1e3


def monthly_means(f, n=12):
    """n equal bins of the (periodic, 1-yr) forcing -> bin means at bin centres"""
    t = np.asarray(f['time'], dtype=float) % 1.0
    b = np.minimum((t * n).astype(int), n - 1)
    out = {'time': (np.arange(n) + 0.5) / n}
    for k, v in f.items():
        if k != 'time':
            v = np.asarray(v, dtype=float)
            out[k] = np.array([v[b == i].mean() for i in range(n)])
    return out


f_month = monthly_means(f_daily)
fields = [k for k in f_month if k != 'time']

cbed2py.set_model(args.model)
print(f'{args.name}: model {args.model}, {len(f_daily["time"])} daily -> '
      f'{len(f_month["time"])} monthly samples, forcing {fields}')
for k in fields:
    print(f'  {k:7s} daily [{f_daily[k].min():.4g}, {f_daily[k].max():.4g}]   '
          f'monthly [{f_month[k].min():.4g}, {f_month[k].max():.4g}]')


# ---------------------------------------------------------------------------
# spin-up with the monthly forcing
# ---------------------------------------------------------------------------
# tol_imbalance is absolute [umol C cm-2 yr-1]; set it relative to J.OM
tol_imb = args.tol_imb * float(np.mean(f_month['J.OM']))
print(f'tol_imbalance = {100*args.tol_imb:g}% of mean J.OM = {tol_imb:.3g} umol C cm-2 yr-1')

t0 = time.time()
y, hist, ds_last = newton_krylov_spinup(
    params, f_month, n_periods=args.n_periods, max_newton=args.max_newton,
    gmres_iter=args.gmres_iter, tol_f=args.tol_f, tol_imbalance=tol_imb,
    checkpoint=f'{args.name}_spinup_monthly_ckpt.npz')
h = hist[-1]
print(f'spin-up {h["runs"] * args.n_periods} model years, {time.time() - t0:.0f} s, '
      f'converged: {h["converged"]}')
if not h['converged']:
    print(f'WARNING: not converged; restart from {args.name}_spinup_monthly_ckpt.npz')

np.savez(f'{args.name}_spinup_monthly.npz', y=y, **{f'forcing_{k}': v for k, v in f_month.items()})


# ---------------------------------------------------------------------------
# check with the daily forcing: the monthly periodic state should be close
# to the daily one, i.e. small OM imbalance and little OC_storage drift
# ---------------------------------------------------------------------------
if args.verify_years > 0:
    n = args.verify_years
    t_exp = np.arange(0, n + 1e-9, 1/365)
    t0 = time.time()
    ds = cbed2py.run_transient(t_exp, forcing=periodic_forcing(f_daily, 0, n),
                               params=params, yini=y)
    print(f'\ncheck run: {n} yr with daily forcing, {time.time() - t0:.0f} s')
    print(' year   imbalance [umol C cm-2 yr-1]   % of J.OM   OC_storage change vs prev. year')
    yrs = np.arange(0, n + 1)
    st = ds.OC_storage.sel(time=yrs, method='nearest').values
    for T in range(1, n + 1):
        imb = annual_budget(ds, T)['imbalance']
        print(f' {T:4d}   {imb:12.4f}                  {100*abs(imb)/params["J.OM"]:7.3f}'
              f'     {(st[T] - st[T-1]) / st[T-1]:+.2e}')

    # -----------------------------------------------------------------------
    # plot: forcing and response over the check run (starts from the spin-up)
    # -----------------------------------------------------------------------
    import matplotlib.pyplot as plt

    FLUX_CONV = 36.5                         # umol cm-2 yr-1 -> mmol m-2 d-1
    C_DAILY, C_MONTH = '#2a78d6', '#eb6834'
    fd = periodic_forcing(f_daily, 0, n, lead=0)
    fm = periodic_forcing(f_month, 0, n, lead=0)

    rows = [('J.OM', 'J.OM forcing', 'umol C cm$^{-2}$ yr$^{-1}$'),
            ('O2.w', 'bottom-water O$_2$', 'mM'),
            ('o2_flux', 'O$_2$ flux', 'mmol m$^{-2}$ d$^{-1}$'),
            ('dic_flux', 'DIC flux', 'mmol m$^{-2}$ d$^{-1}$'),
            ('no3_flux', 'NO$_3$ flux', 'mmol m$^{-2}$ d$^{-1}$'),
            ('nh4_flux', 'NH$_4$ flux', 'mmol m$^{-2}$ d$^{-1}$'),
            ('OC_storage', 'OC storage', 'mol C m$^{-2}$')]

    fig, axs = plt.subplots(len(rows), 1, figsize=(10, 13), sharex=True,
                            constrained_layout=True)
    for ax, (v, title, unit) in zip(axs, rows):
        if v in f_daily:                     # forcing: daily (used) vs monthly (spin-up)
            ax.plot(fd['time'], fd[v], color=C_DAILY, lw=1, label='daily (check run)')
            ax.plot(fm['time'], fm[v], color=C_MONTH, lw=2, label='monthly (spin-up)')
            ax.legend(frameon=False, fontsize=8, loc='upper right')
        else:                                # model response, daily forcing
            y_v = ds[v] / FLUX_CONV if v.endswith('_flux') else ds[v]
            ax.plot(ds.time, y_v, color=C_DAILY, lw=1)
        ax.set_title(title, loc='left', fontsize=10)
        ax.set_ylabel(unit, fontsize=8)
        ax.grid(axis='y', color='0.92', lw=0.6)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
    axs[-1].set_xlabel('time since end of spin-up [yr]')
    fig.suptitle(f'{args.name}: {n} yr with daily forcing from the monthly spin-up state '
                 '(fluxes positive into the sediment)', fontsize=11)
    fig.savefig(f'{args.name}_transient.png', dpi=150)
    print(f'figure: {args.name}_transient.png')


# ---------------------------------------------------------------------------
# drift run: same diagnostics and plots as ../test_krylov_sites.py
# ---------------------------------------------------------------------------
if args.drift_years > 0:
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    FLUX_CONV = 36.5                         # umol cm-2 yr-1 -> mmol m-2 d-1
    n = args.drift_years
    f_drift = f_month if args.drift_forcing == 'monthly' else f_daily
    t_exp = np.arange(0, n + 1e-9, 1/365)
    t0 = time.time()
    ds = cbed2py.run_transient(t_exp, forcing=periodic_forcing(f_drift, 0, n),
                               params=params, yini=y)
    print(f'\ndrift run: {n} yr with {args.drift_forcing} forcing, {time.time() - t0:.0f} s')

    marks = [T for T in (1, 10, 50, 100) if T <= n] + ([n] if n not in (1, 10, 50, 100) else [])
    imb = [annual_budget(ds, T)['imbalance'] for T in marks]
    print(f'annual imbalance at yr {marks} [umol C cm-2 yr-1]:', np.round(imb, 4),
          f'({100*abs(imb[-1])/params["J.OM"]:.3f}% of J.OM)')
    print(f'OC_storage change over {n} yr (relative):',
          float((ds.OC_storage[-1] - ds.OC_storage[0]) / ds.OC_storage[0]))
    tag = f'{args.name} ({args.drift_forcing} forcing)'

    # % change of OM1-3 in the top, mid and bottom layers, sampled once a year
    # at the same phase (t = 0, 1, 2, ...) so the seasonal cycle drops out
    yrs = np.arange(0, n + 1)
    ann = ds.sel(time=yrs, method='nearest')
    N = ds.z.size
    layers = {'top': 0, 'mid': N // 2, 'bottom': N - 1}

    fig, axs = plt.subplots(1, 3, figsize=(11, 3.5), sharex=True, constrained_layout=True)
    for ax, (name, k) in zip(axs, layers.items()):
        for s in ('OM1', 'OM2', 'OM3'):
            c = ann[s].isel(z=k)
            if c[0] > 1e-6:                  # skip pools that are ~0 there (e.g. OM1 at depth)
                ax.plot(yrs, 100 * (c - c[0]) / c[0], label=s)
        ax.axhline(0, color='0.6', lw=0.8)
        ax.set(title=f'{name} layer (z = {float(ds.z[k]):.2f} cm)', xlabel='time [yr]')
        ax.legend(frameon=False)
    axs[0].set_ylabel('change vs year 0 [%]')
    fig.suptitle(f'{tag}: drift after spin-up')
    fig.savefig(f'{args.name}_drift_layers.png', dpi=150)

    # Hovmoller (time x depth) of the % change of OM1-3 vs year 0, same annual
    # sampling. Each panel has its own symmetric colour range; layers where
    # the pool is ~0 are masked.
    fig, axs = plt.subplots(1, 3, figsize=(12, 4), sharey=True, constrained_layout=True)
    for ax, s in zip(axs, ('OM1', 'OM2', 'OM3')):
        c = ann[s]                                # (time, z)
        c0 = c.isel(time=0)
        pct = 100 * (c - c0) / c0.where(c0 > 1e-3 * c0.max())
        vmax = float(np.nanmax(np.abs(pct))) or 1e-12
        pm = ax.pcolormesh(yrs, ds.z, pct.T, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                           shading='nearest')
        fig.colorbar(pm, ax=ax, label='change vs year 0 [%]')
        ax.set(title=s, xlabel='time [yr]')
    axs[0].invert_yaxis()                         # shared y: inverts all panels
    axs[0].set_ylabel('depth [cm]')
    fig.suptitle(f'{tag}: OM drift after spin-up (same phase each year)')
    fig.savefig(f'{args.name}_drift_hovmoller.png', dpi=150)

    # Hovmoller (time x depth) of the OM1-3 concentrations, full output
    # resolution (seasonal cycle included). Log colour scale; values ~0 masked.
    fig, axs = plt.subplots(3, 1, figsize=(12, 6), sharey=True, constrained_layout=True)
    for ax, s in zip(axs, ('OM1', 'OM2', 'OM3')):
        c = ds[s]                                 # (time, z)
        cmax = float(c.max())
        cm = c.where(c > 1e-6 * cmax)
        pm = ax.pcolormesh(ds.time, ds.z, cm.T, cmap='viridis', shading='nearest',
                           norm=LogNorm(vmin=float(cm.min()), vmax=cmax))
        fig.colorbar(pm, ax=ax, label=f'{s} [umol cm-3]')
        ax.set(title=s, xlabel='time [yr]')
    axs[0].invert_yaxis()
    axs[0].set_ylabel('depth [cm]')
    fig.suptitle(f'{tag}: OM concentrations')
    fig.savefig(f'{args.name}_om_hovmoller.png', dpi=150)

    # benthic fluxes, umol cm-2 yr-1 -> mmol m-2 d-1; positive into the sediment
    fig, axs = plt.subplots(4, 1, figsize=(12, 6), sharex=True, constrained_layout=True)
    for ax, s in zip(axs, ('o2_flux', 'dic_flux', 'nh4_flux', 'no3_flux')):
        ax.plot(ds.time, ds[s] / FLUX_CONV)
        ax.axhline(0, color='0.6', lw=0.8)
        ax.set(title=s, ylabel='mmol m$^{-2}$ d$^{-1}$')
    axs[-1].set_xlabel('time [yr]')
    fig.suptitle(f'{tag}: benthic fluxes (positive into the sediment)')
    fig.savefig(f'{args.name}_drift_fluxes.png', dpi=150)
    print(f'figures: {args.name}_drift_layers.png, _drift_hovmoller.png, '
          f'_om_hovmoller.png, _drift_fluxes.png')

if args.verify_years > 0 or args.drift_years > 0:
    import matplotlib.pyplot as plt
    plt.show()
