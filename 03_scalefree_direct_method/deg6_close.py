import sympy as sp, numpy as np, sys, time
from sf_general import build, th, S, C, redSC, J

def setup(D=10):
    params = sp.symbols('a0:3')
    eqsJ, Uj, unknowns = build(6, 2, params, 1)
    # trig-polynomial ansatz for each unknown function
    coeffs = []; rep = {}
    for key in Uj:
        expr = 0
        for i in range(D+1):
            for j in range(2):
                if i + j > D: continue
                sym = sp.Symbol(f'u{key[0]}_{key[1]}_{i}_{j}'); coeffs.append(sym)
                expr += sym*sp.sin(th)**i*sp.cos(th)**j
        for k in range(J-1, -1, -1):
            rep[Uj[key][k]] = sp.diff(expr, th, k)
    return eqsJ, params, coeffs, rep

def system_at(eqsJ, params, coeffs, rep, c2, c4, K=1):
    c0 = 1
    # Killing params from (c0,c2,c4): c0 = (45a0+15a1+45a2)/120, c2 = 20(a0-a2)/120, c4 = (a1-a0-a2)/120  (K=1)
    a0, a1, a2 = params
    sol = sp.solve([sp.Rational(1,120)*(45*a0+15*a1+45*a2) - c0, sp.Rational(1,6)*(a0-a2) - c2, sp.Rational(1,120)*(a1-a0-a2) - c4], params)
    fth = c0 + c2*sp.cos(2*th) + c4*sp.cos(4*th)
    Fsub = {sp.Symbol(f'F{k}'): sp.expand(sp.expand_trig(sp.diff(fth, th, k))) for k in range(J)}
    lin = []
    for key, e in eqsJ.items():
        ex = e.subs(sol).subs(Fsub).subs(rep)
        ex = sp.expand(ex.subs({sp.sin(th): S, sp.cos(th): C}))
        ex = redSC(ex)
        P = sp.Poly(ex, S, C)
        lin += list(P.coeffs())
    M, rhs = sp.linear_eq_to_matrix(lin, coeffs)
    return M, rhs

if __name__ == '__main__':
    t0 = time.time()
    eqsJ, params, coeffs, rep = setup()
    print("setup", time.time()-t0, "unknown coefficients:", len(coeffs), flush=True)
    # exact checks at known points
    for c2, c4, label in [(0, 0, 'sphere'), (sp.Rational(1,3), 0, 'parabolic (z)'), (-sp.Rational(1,3), 0, 'parabolic (x)'), (sp.Rational(3,10), 0, 't=0.3 control')]:
        M, rhs = system_at(eqsJ, params, coeffs, rep, c2, c4)
        Mn = np.array(M, dtype=float); rn = np.array(rhs, dtype=float).ravel()
        x, res, rk, sv = np.linalg.lstsq(Mn, rn, rcond=None)
        r = np.linalg.norm(Mn @ x - rn)
        print(f"  {label:18s} c2={float(c2):+.3f} c4={float(c4):+.3f}: residual {r:.2e}", flush=True)
    # scan
    grid = np.linspace(-1, 1, 41)
    res = np.zeros((41, 41))
    for i, c2 in enumerate(grid):
        for j, c4 in enumerate(grid):
            M, rhs = system_at(eqsJ, params, coeffs, rep, sp.nsimplify(round(c2, 3)), sp.nsimplify(round(c4, 3)))
            Mn = np.array(M, dtype=float); rn = np.array(rhs, dtype=float).ravel()
            x, _, _, _ = np.linalg.lstsq(Mn, rn, rcond=None)
            res[i, j] = np.linalg.norm(Mn @ x - rn)/max(1.0, np.linalg.norm(rn))
        print(f"  row {i} done ({time.time()-t0:.0f}s), min so far {res[:i+1].min():.2e}", flush=True)
    np.save('deg6_d2_scan.npy', res)
    hits = [(grid[i], grid[j], res[i, j]) for i in range(41) for j in range(41) if res[i, j] < 1e-8]
    print("hits:", hits)
    from scipy.ndimage import minimum_filter
    loc = (res == minimum_filter(res, size=5)) & (res < 1e-2)
    print("local minima < 1e-2:", [(round(grid[i], 3), round(grid[j], 3), '%.1e' % res[i, j]) for i, j in zip(*np.where(loc))])
