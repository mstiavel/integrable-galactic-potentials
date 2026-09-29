import sys, sympy as sp, numpy as np, pickle
from scipy.optimize import least_squares
from second_order import run
N, M = int(sys.argv[1]), int(sys.argv[2])
conds, tpar, lowb, sub1, b1 = run(N, M, 1, 0, verbose=True)
# normalise: b2 = 1  (lowb[1] = t9/6 for beta=0)
tsyms = list(tpar)
f = sp.lambdify(tsyms, conds, 'numpy')
b2expr = sp.lambdify(tsyms, [lowb[1]], 'numpy')
def resid(t):
    return np.array(f(*t), dtype=float).tolist() + [10.0*(b2expr(*t)[0] - 1.0)]
rng = np.random.default_rng(0)
best = None
for trial in range(60):
    t0 = rng.normal(size=len(tsyms))*2
    r = least_squares(resid, t0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
    if best is None or r.cost < best.cost: best = r
    if r.cost < 1e-24: break
print(f"N={N} M={M}: best residual^2 = {best.cost:.3e}  (with b2 fixed to 1)")
t = best.x
vals = [float(sp.lambdify(tsyms, x, 'numpy')(*t)) for x in lowb]
print("   B1 coefficients b1..:", np.round(vals, 8))
pickle.dump((best.x, tsyms), open(f'so_num_{N}_{M}.pkl','wb'))
