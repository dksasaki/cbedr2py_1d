# Realistic forcing: spin-up for one grid point

Newton-Krylov spin-up of `cbed2py` to its periodic state for a single model
grid point, with bottom-water forcing taken from a model. The forcing is not
synthetic: it is one year of daily output from that model.

## Files

| file | content |
|---|---|
| `<name>_params.json` | params dict (CBED units), written by `extract_forcing.ipynb` |
| `<name>_forcing.npz` / `.csv` | one year of daily forcing: `time` [yr] + `J.OM`, `O2.w`, `NO3.w`, `DIC.w`, `NH4.w`, `TAlk.w` |
| `spinup_monthly.py` | **recommended**: spin-up with monthly-mean forcing, then a check run with the daily forcing |
| `spinup_from_forcing.py` | spin-up with the daily forcing as is (slow, see below) |

The example point is `forcing_ix294_iy100`.

## Run

From this folder, inside the pixi environment:

```bash
pixi run python spinup_monthly.py forcing_ix294_iy100
pixi run python spinup_monthly.py forcing_ix294_iy100 --tol-imb 0.005     # stricter budget
pixi run python spinup_monthly.py forcing_ix294_iy100 --verify-years 0    # skip check run
```

| option | default | meaning |
|---|---|---|
| `--model` | `fast` | `fast` or `v2` |
| `--n-periods` | 5 | years per Newton shot (use 10 for slow, low-flux sites) |
| `--max-newton` | 10 | maximum Newton steps |
| `--gmres-iter` | 20 | maximum GMRES iterations per Newton step |
| `--tol-f` | 3e-5 | relative change of the state per period (RMS) |
| `--tol-imb` | 0.02 | OM budget imbalance, as a fraction of mean J.OM |
| `--verify-years` | 3 | years of the check run with daily forcing |

Outputs:
- `<name>_spinup_monthly.npz`: the periodic state `y` (cbed2py state order)
  and the monthly forcing used to reach it.
- `<name>_spinup_monthly_ckpt.npz`: a checkpoint written after every Newton
  step. If the run does not converge, restart from it.

## Why monthly forcing

The R model interpolates the forcing linearly in time. Every forcing sample
is therefore a kink in the boundary conditions, and lsoda shortens its time
step at each one. Smooth forcing (a sine, as in `../test_krylov_sites.py`)
barely notices the kinks. Noisy daily bottom-water data costs much more,
especially `O2.w`, because O2 drives the reactions.

Cost of one model year (fast model, this point):

| forcing | s / model yr |
|---|---|
| daily, all fields | 7.6 |
| daily, only `O2.w` | 7.1 |
| daily, only `J.OM` | 2.8 |
| weekly samples | 2.2 |
| monthly samples | 0.6 |
| smooth sine `J.OM` / constant | 0.2 |

A spin-up needs about 100 model years or more, so daily forcing takes tens of
minutes and monthly forcing takes about a minute. Daily noise in the bottom
water has little effect on the slow OM pools that set the periodic state.
The check run measures how big that effect is.

## Convergence criteria

The spin-up stops when both criteria are met:

- `tol_f`: the RMS over all cells and species of (change per period / state
  size). This is relative: 3e-5 means about 0.003 % change per period. Do
  not go below about 1e-5, because that is the noise floor of lsoda.
- `tol_imbalance`: the annual-mean OM_flux − remineralisation − burial.
  `newton_krylov_spinup` takes it in absolute units [µmol C cm⁻² yr⁻¹].
  `spinup_monthly.py` sets it as `--tol-imb` × mean J.OM, so that the same
  setting means the same relative accuracy at every site.

## Result for `forcing_ix294_iy100`

Settings: monthly forcing, defaults, with `tol_imbalance` = 0.05 absolute
(before `--tol-imb` was added).

- **Spin-up:** 5 Newton steps, 120 model years, 61 s. Final imbalance 0.034
  (0.04 % of J.OM).
- **Check run (3 years, daily forcing):**
  - Imbalance about 0.15 µmol C cm⁻² yr⁻¹ (0.17 % of J.OM).
  - OC_storage change about +2e-5 per year.

The periodic state for monthly forcing is within about 0.2 % of the one for
daily forcing. If that is not enough, pass the state as `yini` to
`newton_krylov_spinup` with the daily forcing. It should then need only one
or two Newton steps.

## Known issue: TAlk units

`TAlk.w` comes out of `extract_forcing.ipynb` in µmol/kg (about 2414), but the
model expects mM (about 2.41). `spinup_monthly.py` detects values above 100
and divides them by 1000, and prints a note when it does. `spinup_from_forcing.py`
does not. TAlk is passive, so it affects neither the run time nor the other
species, but the TAlk profiles and `talk_flux` are wrong by a factor of 1000.
The proper fix belongs in the notebook.
