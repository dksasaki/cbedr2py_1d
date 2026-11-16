"""
Newton-Krylov spin-up of cbed2py from a saved point forcing (no project code needed).

Input: the files written by extract_forcing.ipynb
    <name>_params.json   params dict
    <name>_forcing.npz   f_year dict (time + J.OM, O2.w, ...)

Usage (from this folder):
    python spinup_from_forcing.py forcing_ix294_iy100
    python spinup_from_forcing.py forcing_ix294_iy100 --n-periods 2 --max-newton 6
    python spinup_from_forcing.py forcing_ix294_iy100 --model v2

Same call as the site script (newton_krylov_spinup(params, f_year, ...)),
only params and f_year come from the files. Writes <name>_spinup.npz
(the periodic state y, in cbed2py state order) and a checkpoint.
"""
import argparse
import json
import time
import numpy as np
import cbed2py
from cbed2py.spinup import newton_krylov_spinup

ap = argparse.ArgumentParser()
ap.add_argument('name', help='file stem, e.g. forcing_ix294_iy100')
ap.add_argument('--model', default='fast', help="'fast' or 'v2'")
ap.add_argument('--n-periods', type=int, default=5)
ap.add_argument('--max-newton', type=int, default=10)
ap.add_argument('--gmres-iter', type=int, default=20)
ap.add_argument('--tol-f', type=float, default=3e-5)
args = ap.parse_args()

with open(f'{args.name}_params.json') as f:
    params = json.load(f)
d = np.load(f'{args.name}_forcing.npz')
f_year = {k: d[k] for k in d.files}          # 'time' + forcing fields

cbed2py.set_model(args.model)
print(f'{args.name}: model {args.model}, {len(f_year["time"])} phases, '
      f'forcing {[k for k in f_year if k != "time"]}')

t0 = time.time()
y, hist, ds_last = newton_krylov_spinup(
    params, f_year, n_periods=args.n_periods, max_newton=args.max_newton,
    gmres_iter=args.gmres_iter, tol_f=args.tol_f,
    checkpoint=f'{args.name}_spinup_ckpt.npz')
h = hist[-1]
print(f'spin-up {h["runs"] * args.n_periods} model years, {time.time() - t0:.0f} s, '
      f'converged: {h["converged"]}')
if not h['converged']:
    print(f'WARNING: not converged; restart from the checkpoint {args.name}_spinup_ckpt.npz')

np.savez(f'{args.name}_spinup.npz', y=y)
