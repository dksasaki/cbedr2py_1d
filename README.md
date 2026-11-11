# cbed2py_1d — transient CBED at a single point

`cbed2py` is a Python package that wraps (with rpy2) a **time-evolving**
version of the CBED v2 sediment model. It has the same equations and
coefficients as `cbed_model()` in `cbed_R2py/src/cbed_v2_func.R`, the
steady-state model used by `cbed_R2py/scripts/cbed_wrapper.py`. The forcing
can now vary in time. The project is independent of `cbed_R2py`: it has its
own pixi environment, and the R files it needs are copied into
`cbed2py/model/`.

## Contents

1. [Files](#1-files)
2. [Installing and running](#2-installing-and-running)
3. [The wrapper in one example](#3-the-wrapper-in-one-example)
4. [Functions](#4-functions)
5. [Inputs: what goes where, and in which format](#5-inputs-what-goes-where-and-in-which-format)
6. [Parameters](#6-parameters)
7. [Output](#7-output)
8. [From MOM6-COBALT to CBED units](#8-from-mom6-cobalt-to-cbed-units)
9. [Differences from the steady-state model](#9-differences-from-the-steady-state-model)
10. [Checks](#10-checks)

---

## 1. Files

```
cbed2py_1d/
├── pixi.toml                     environment (conda-forge + CRAN) and tasks
├── pyproject.toml                makes cbed2py installable; pixi installs it in editable mode
├── README.md
├── cbed2py/                      the package
│   ├── __init__.py               exposes default_params, steady_state, run_transient, ...
│   ├── wrapper.py                the wrapper
│   └── model/
│       ├── cbed_v2_transient.R   the transient model
│       ├── default_CBED_pars.R   get_default_pars()
│       └── cbed_v2_func.R        the original steady-state cbed_model() (reference)
└── tests/
    └── test_wrapper.py           four checks (section 10)
```

| file | role |
|---|---|
| `cbed2py/wrapper.py` | the wrapper: `default_params()`, `steady_state()`, `run_transient()` |
| `cbed2py/model/cbed_v2_transient.R` | the transient model: `cbed_model_t()`, `cbed_model_ss_t()`, `cbed_t_default_coefs()` |
| `cbed2py/model/default_CBED_pars.R` | `get_default_pars()`, copied from `cbed_R2py/scripts/` |
| `cbed2py/model/cbed_v2_func.R` | the original steady-state `cbed_model()`; only `tests/test_wrapper.py` uses it, as a reference |
| `tests/test_wrapper.py` | four checks (section 10) |

## 2. Installing and running

Requires [pixi](https://pixi.sh). All dependencies are in `pixi.toml`, with
the same versions as `cbed_R2py`. `cbed2py` is installed in editable mode, so
changes to the code take effect without reinstalling.

```bash
cd cbed2py_1d
pixi install          # environment in .pixi/, installs cbed2py
pixi run setup-r      # once: ReacTran and marelac from CRAN (not on conda-forge)
pixi run test         # the five checks
```

Then use it with `pixi run python my_script.py`, or open a shell in the
environment with `pixi shell`.

R starts the first time a wrapper function is called. A steady state takes
about 30 s. A 2-year run with daily output takes about 25 s. Memory use is
about 300 MB.

## 3. The wrapper in one example

```python
import numpy as np
import cbed2py

# point values and constant (spin-up) forcing, CBED units
params = {
    'por.0': 0.75, 'por.inf': 0.75,          # porosity
    'depth': 250.0,                          # m
    'f.OM1': 0.33, 'f.OM2': 0.62, 'f.OM3': 0.05,
    'J.OM': 400.0, 'O2.w': 0.25, 'NO3.w': 0.015, 'NH4.w': 0.001,
    'DIC.w': 2.15, 'TAlk.w': 2.25, 'temp': 7.0, 'S': 33.0, 'w': 0.05,
}

# time-varying forcing, CBED units, time in years
t = np.arange(0, 3 + 1e-9, 1/365)
forcing = {'time': t,
           'J.OM': 400 * (1 + 0.5*np.sin(2*np.pi*t)),
           'temp': 7 + 2*np.sin(2*np.pi*(t - 0.6))}

out = cbed2py.run_transient(t, forcing=forcing, params=params)   # spin-up from steady state
out.to_netcdf('example.nc')
```

## 4. Functions

### `default_params() -> dict`
Returns every parameter with its default value: `get_default_pars()` plus the
v2 coefficients from `cbed_t_default_coefs()` (section 6). Use it to see
what can be changed.

### `steady_state(params=None) -> xr.Dataset`
Computes the steady state of the transient equations, using the constant
forcing given in `params`. `params` is a dict of overrides applied on top of
`default_params()`; unknown names raise a `KeyError`. Returns the 9 state
profiles on `z` and the flux diagnostics as scalars. `attrs['steady']` is
`1` if the solver converged.

### `run_transient(times, forcing=None, params=None, yini=None) -> xr.Dataset`

| argument | what | format |
|---|---|---|
| `times` | output times [yr]; the first value is the initial time | 1-D array |
| `forcing` | time series of boundary conditions | `dict` or `pandas.DataFrame` (section 5) |
| `params` | overrides of `default_params()`: point values, spin-up forcing, coefficients | `dict` of floats |
| `yini` | initial state | `None`, `xr.Dataset` or array of 9·N |

`yini`:
- `None` (the default): `steady_state(params)` is computed first. The
  constant forcing values in `params` (`J.OM`, `O2.w`, …) set the spin-up.
  Use the time mean of the forcing for them.
- A Dataset from `steady_state()` or `run_transient()`. For a run, its last
  time is used, so passing the previous output continues that run.
- An array of `9*N` values, in the order OM1…TAlk, each block running from
  the surface down.

If the solver stops early, a warning is printed and the output has fewer
times than requested.

## 5. Inputs: what goes where, and in which format

CBED needs two kinds of input, and they go in different arguments:

| kind | examples | can vary in time? | goes in |
|---|---|---|---|
| **boundary forcing** | `J.OM`, `O2.w`, `NO3.w`, `NH4.w`, `DIC.w`, `TAlk.w`, `temp`, `S`, `w` | yes | `forcing` (time series) **and** `params` (the spin-up value) |
| **point properties** | `por.0`, `por.inf`, `por.att`, `depth`, `f.OM1-3` | no | `params` (scalars) |
| **coefficients** | `k1.0`, `k2.0`, `Db.coef`, … | no | `params` (scalars) |

**Porosity is a point property, not a forcing.** The porosity profile, the
solid volume fraction and the grid are built once, before the integration.
They can't change during a run. Set `por.0` (surface) and `por.inf` (at
depth) in `params`; the steady-state wrapper sets both to the same value.
`por.att` sets the attenuation length (cm) between the two. Porosity also
enters the forcing through `w`: the caller computes `w = w_tem / (1 - por)`
(section 8), so use the same porosity in both places.

### Forcing format today

`forcing` is passed through `pandas.DataFrame(forcing)`, so the following
work:

| format | works? | how |
|---|---|---|
| `dict` of 1-D arrays | yes | `{'time': t, 'J.OM': a, ...}` |
| `pandas.DataFrame` | yes | columns `time`, `J.OM`, … |
| `xr.Dataset` | via pandas | `ds.to_dataframe().reset_index()`, with `time` as float years |
| `xr.DataArray` | no | one per variable; put them in a dict or Dataset first |

Rules (checked by the wrapper):
- a `time` column is required. It must be **float years**, strictly
  increasing, and on the same time axis as `times`.
- other columns must be among the 9 forcing names above. Anything else
  raises a `KeyError`, which catches typos such as `O2_w`.
- columns can be left out: those variables stay constant at their `params`
  value.
- values must already be in **CBED units** (section 8). The wrapper converts
  nothing.
- between samples the model interpolates linearly. Before the first and after
  the last sample the value is held constant (R `approxfun(rule = 2)`).
  Forcing does not repeat on its own: to cycle a year, tile it with shifted
  times.

### Recommendation

- **Forcing → `xr.Dataset`.** The MOM6-COBALT output is read with xarray
  anyway. A Dataset with a `time` coordinate and one variable per CBED
  forcing (named `J.OM`, `O2.w`, …) is the natural object to build, check,
  plot and save to netCDF before running. Units and the source of each
  variable can go in `attrs`. Today it has to be flattened with
  `to_dataframe().reset_index()`, and time has to be converted to years by
  hand. A small change to `run_transient()` would let it accept the Dataset
  directly, converting a `datetime64` time to years since the first record.
  That change isn't made yet.
- **Porosity and the other point properties → scalars in `params`.** They
  are single numbers per point, and the model can't use a time series for
  them. If they come from a file, read the value at the point (`float(...)`)
  and put it in the dict.
- **Coefficients → `params`.** Only list the ones you change; the defaults
  fill in the rest.
- Plain `dict`s of arrays are fine for synthetic or quick tests.

## 6. Parameters

`default_params()` merges two lists.

**From `get_default_pars()` (`cbed2py/model/default_CBED_pars.R`)**, the ones you will
normally set:

| name | meaning | unit | default |
|---|---|---|---|
| `J.OM` | POC deposition flux | µmol C cm⁻² yr⁻¹ | 403.3 |
| `O2.w`, `NO3.w`, `NH4.w` | bottom-water concentrations | mM (= mol m⁻³) | 0.276, 0.008, 0.0006 |
| `DIC.w`, `TAlk.w` | bottom-water DIC, alkalinity | mM | 2.16, 2.2 |
| `ODU.w` | bottom-water ODU (not forced) | mM | 0 |
| `temp`, `S` | bottom temperature, salinity | °C, psu | 6.87, 34 |
| `w` | sedimentation rate at the surface | cm yr⁻¹ | 0.0359 |
| `por.0`, `por.inf`, `por.att` | porosity at the surface and at depth; attenuation length | –, –, cm | 0.67, 0.67, 2 |
| `depth` | water depth (`f.OM` only) | m | 74 |
| `f.OM1`, `f.OM2`, `f.OM3` | fractions of `J.OM` into the 3 OM pools; must add to 1 | – | 0.7, 0.2, 0.1 |
| `N`, `L` | number of layers, sediment depth | –, cm | 20, 20 |

Also in that list: half-saturations (`Ks.O2`, `Ks.NO3`), reaction rates
(`k.NOx`, `k.ana`, `k.ODUox`), `Q10`, `K.equilNH4`, `ps`, `rNC*`,
`irr.enh.*` and `switch.odu.burial`. Some entries are unused (`Db.max`,
`alpha.max`, `xL`, `switch.oxy.dependent.k`, `run_*`).

`depth` doesn't enter the equations. The steady-state wrapper uses it to set
the pool fractions:
`f.OM1 = min(0.80, 0.65·(100/depth)^0.5)`,
`f.OM3 = max(0.03, 0.04·(100/depth)^-0.3)`, `f.OM2 = 1 − f.OM1 − f.OM3`.
The transient wrapper does not compute them, so pass `f.OM1-3` yourself.

**From `cbed_t_default_coefs()` (`cbed2py/model/cbed_v2_transient.R`):** constants
that are hard-coded in `cbed_v2_func.R`, with the v2 values as defaults.

The three `jpoc_mode_*` parameters pick the POC flux `J` used by reactivity,
bioturbation and bioirrigation, separately. Each is one of the strings
`'instantaneous'` (the forcing `J.OM(t)`, default), `'longterm_mean'` (the
constant `params['J.OM']`, the role of the Fortran `cbed_jpoc_mean`) or
`'trailing_mean'` (the mean of `J.OM(t)` over the previous `J.avg.days`
days). Any other value raises a `ValueError`. The wrapper passes them to R as
integer codes 0/1/2 under the R names `J.mode.k`, `J.mode.Db`, `J.mode.irr`.
The Fortran flags `use_jpoc_mean_for_k = .true.` /
`use_jpoc_mean_for_mixing = .true.` (its defaults) correspond to
`'longterm_mean'`.

| name | default | used in |
|---|---|---|
| `k1.0`, `k2.0`, `k3.0` | 0.15, 1.5e-3, 0.9e-4 | `k_i = k_i.0 · J^k.exp` [yr⁻¹] |
| `k.exp` | 0.85 | same |
| `k.adj.denit`, `k.adj.anoxia` | 0.1, 0.001 | decay factor under denitrification / anoxia |
| `jpoc_mode_reactivity` | `'instantaneous'` | `J` in `k_i` |
| `jpoc_mode_bioturbation` | `'instantaneous'` | `J` in `Db0` |
| `jpoc_mode_bioirrigation` | `'instantaneous'` | `J` in `I0` |
| `J.avg.days` | 7 d | window of `'trailing_mean'` |
| `Db.coef`, `Db.exp` | 0.0232, 0.85 | `Db0 = Db.coef · J^Db.exp` |
| `Db.xL`, `Db.O2.half` | 8 cm, 20 µM | `Db = Db0 · exp(−(z/Db.xL)²) · O2/(O2 + Db.O2.half)` |
| `irr.xL` | 1.8 cm | `alpha = I0 · exp(−(z/irr.xL)²)` |
| `odu.bur.a`, `odu.bur.b` | 0.233, 0.336 | ODU burial `min(1, a · w^b)` |

`k1.0`–`k3.0` and `k.adj.*` match `generic_CBED.F90` (lines 1134–1135 and
1324–1326).

## 7. Output

`run_transient()` returns an `xr.Dataset`:

| variable | dims | unit |
|---|---|---|
| `OM1`, `OM2`, `OM3` | (time, z) | µmol cm⁻³ of solid |
| `O2`, `NH4`, `NO3`, `ODU`, `DIC`, `TAlk` | (time, z) | mM |
| `o2_flux`, `nh4_flux`, `no3_flux`, `dic_flux`, `odu_flux`, `talk_flux` | (time) | µmol cm⁻² yr⁻¹, positive into the sediment |
| `OM_flux`, `OM_burial`, `tot_remin`, `oxic_remin`, `tot_denit`, `anammox`, `nitrification` | (time) | µmol cm⁻² yr⁻¹ |
| `OC_storage` | (time) | as in `cbed_v1_func.R` |
| `Db_top`, `alpha_top` | (time) | bioturbation, bioirrigation at the top layer |
| `forcing_<name>` | (time) | the forcing the model saw, at the output times |

- Coordinates: `time` [yr], `z` [cm, layer mid-points].
- All parameters used are stored as attributes (`par_<name>`).
- `steady_state()` returns the same states on `z` only, with the fluxes as
  scalars.
- To convert fluxes to mmol m⁻² d⁻¹, multiply by `10/365`.

## 8. From MOM6-COBALT to CBED units

These are the formulas `cbed_R2py/scripts/cbed_wrapper.py` (`cbed_wrapped()`)
uses. The caller applies them before building the forcing. `max(0, ·)` is
applied to each result.

| CBED | from COBALT | conversion |
|---|---|---|
| `J.OM` | `fntot = fndet_btm + fndi_btm + fnsm_btm + fnmd_btm + fnlg_btm` [mol N m⁻² s⁻¹] | `fntot · 6.625 · 1000 · 86400 · 365 / 10` |
| `O2.w` | `btm_o2` [mol kg⁻¹] | `· 1e3` |
| `NO3.w` | `btm_no3` [mol kg⁻¹] | `· 1e3` |
| `NH4.w` | `nh4`, deepest valid level of `cobalt_tracers` [mol kg⁻¹] | `· 1e3` |
| `DIC.w` | `btm_dic` [mol kg⁻¹] | `· 1e3` |
| `TAlk.w` | `talk`, deepest valid level [mol m⁻³] | none |
| `temp`, `S` | `temp`, `salt`, bottom level of `ocean_daily` | none |
| `w` | particle fluxes (below) | `w_tem / (1 − porosity)` |

```
w_tem [cm/yr] = ( fcadet_arag_btm·100/2.71 + fcadet_calc_btm·100/2.94
                + (fsidet_btm + fsimd_btm + fsilg_btm)·60/2.65
                + (ffedet_btm + ffedi_btm + ffesm_btm + ffemd_btm + ffelg_btm)·160/5.24
                + (fpdet_btm + fpdi_btm + fpsm_btm + fpmd_btm + fplg_btm)·120/2.3
                + flithdet_btm/2.65
                + fntot·6.625·22.4/0.9 ) / 1e4 · 86400·365
```

The COBALT `fntot_btm` in `generic_COBALT.F90` (line 3028) also includes
fast-sinking detritus (`f_ndet_fast_btf`). Fast-sinking detritus is not
active in the current COBALT setup, so leaving it out of the sum above is
correct. Add it if it is switched on.

In this output set the files have different sampling: `cobalt_btm` is daily,
`cobalt_tracers` 4-daily, and `ocean_daily` monthly means. Put them on one
time axis before building `forcing`.

## 9. Differences from the steady-state model

`cbed2py/model/cbed_v2_transient.R` started from `cbed_model.t()` in
`cbed_R/cbed_v1_func.R` and was changed to the v2 equations:

- v2 rate constants: `k2 = 1.5e-3`, `k3 = 0.9e-4`, `k.adj.anoxia = 0.001`.
  These are now parameters.
- `f.OM1-3` come from `params`.
- `R.TAlk` uses `−1·R.ODUox`.
- ODU burial is scaled by svf/por.
- Compaction uses `v.0 = w`.
- Forcing comes from the `forcing` argument instead of global `bc.data`.
- `switch.oxy.dependent.k` is ignored: the rate constants always come from
  `k1.0`–`k3.0` and `k.adj.*` (v2's active branch, which is also the Fortran
  formulation). v2's inactive branch (`k3 = 0`) is not reproduced.
- **POC flux for reactivity and mixing.** By default `k1`–`k3`,
  bioturbation and bioirrigation all follow `J.OM(t)`. Each can instead use
  the long-term mean (`'longterm_mean'`, the constant `params['J.OM']`) or
  the trailing mean over the previous `J.avg.days` days (`'trailing_mean'`),
  through `jpoc_mode_reactivity`, `jpoc_mode_bioturbation` and
  `jpoc_mode_bioirrigation`. The trailing mean is the exact mean of the
  linearly interpolated forcing. Before the first forcing time, `J.OM` is
  held at its first value, as for the instantaneous forcing. With constant
  forcing (and in the steady state) all three modes give `params['J.OM']`.
  The `J` actually used is returned as the diagnostics `J_k`, `J_Db`,
  `J_irr`. The deposition flux always follows `J.OM(t)`.

With constant forcing, the transient equations have the same steady state as
`cbed_model()`.

## 10. Checks

`pixi run test` (about 1.5 min). Results on 2026-11-09, default
parameters:

| check | result |
|---|---|
| 1. steady state vs `cbed_model()` of `cbed_v2_func.R` | max relative difference ≤ 8.5e-10 on all 9 profiles; `o2_flux` 379.542, `dic_flux` −342.368, `OM_burial` 60.861 in both |
| 2. constant forcing for 1 yr, from the steady state | max relative drift ≤ 2.9e-11 |
| 3. `J.OM` · (1 + 0.5 sin 2πt), 2 yr daily, default modes (`'instantaneous'`) | 26 s; `OM_flux` follows the forcing (202–605); `o2_flux` 268–474; `dic_flux` −440 to −234; `tot_denit` 16.4–20.5 |
| 4. `k2.0` × 2 | surface OM2 2222 → 1360; `o2_flux` 379.5 → 389.3 |
| 5. `J` modes, monthly `J.OM` cycle, 3 months daily: reactivity `'trailing_mean'` (7 d), bioturbation `'longterm_mean'`, bioirrigation `'instantaneous'` | `J_k` matches a numpy 7-day average to 7e-5 (reference grid resolution); `J_Db` = `J.OM`, `J_irr` = `J.OM(t)` exactly; `J.OM(t)` 202–605, trailing mean 220–587 |

## 11. Disclaimer
Parts of the wrapper were written with assistance from Claude (Anthropic).