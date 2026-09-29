import sympy as sp, numpy as np, time
from deg6_close import setup, system_at
from scipy.ndimage import minimum_filter
t0 = time.time()
eqsJ, params, coeffs, rep = setup()
c2, c4 = sp.symbols('c2 c4')
M, rhs = system_at(eqsJ, params, coeffs, rep, c2, c4)
print("symbolic system built", time.time()-t0, M.shape, flush=True)
Mf = sp.lambdify((c2, c4), M, 'numpy'); rf = sp.lambdify((c2, c4), rhs, 'numpy')
def resid(x2, x4):
    Mn = np.array(Mf(x2, x4), dtype=float); rn = np.array(rf(x2, x4), dtype=float).ravel()
    x, _, _, _ = np.linalg.lstsq(Mn, rn, rcond=None)
    return np.linalg.norm(Mn @ x - rn)/max(1.0, np.linalg.norm(rn))
for x2, x4, lab in [(0,0,'sphere'),(1/3,0,'parabolic z'),(-1/3,0,'parabolic x'),(0.3,0,'t=0.3')]:
    print(f"  {lab:12s}: {resid(x2,x4):.2e}", flush=True)
grid = np.linspace(-1, 1, 81)
res = np.array([[resid(a, b) for b in grid] for a in grid])
np.save('deg6_d2_scan.npy', res)
print("scan done", time.time()-t0, "min residual", res.min(), flush=True)
loc = (res == minimum_filter(res, size=5)) & (res < 1e-2)
print("local minima < 1e-2:", [(round(grid[i],3), round(grid[j],3), '%.1e' % res[i,j]) for i, j in zip(*np.where(loc))])
hits = [(grid[i], grid[j]) for i in range(81) for j in range(81) if res[i,j] < 1e-8]
print("exact hits:", hits)
