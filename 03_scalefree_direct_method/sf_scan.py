import sympy as sp, numpy as np, sys, itertools
from sf_direct import system, th, f
from sf_linalg import canon, S, C
t2, t4 = sp.symbols('t2 t4')

def build(w, D=8, muval=1):
    nparams = {2: 3, 4: 2, 6: 1}[w]
    params = sp.symbols(f'a0:{nparams}')
    eqs, (b20, b11, b02, c) = system(w, params, muval)
    P2 = (3*sp.cos(th)**2 - 1)/2; P4 = (35*sp.cos(th)**4 - 30*sp.cos(th)**2 + 3)/8
    fexpr = 1 + t2*P2 + t4*P4
    coeffs = []
    def ansatz(name):
        terms = 0
        for i in range(D+1):
            for j in range(2):
                if i + j > D: continue
                sym = sp.Symbol(f'{name}_{i}_{j}'); coeffs.append(sym)
                terms += sym*sp.sin(th)**i*sp.cos(th)**j
        return terms
    A20, A11, A02, Ac = ansatz('b20'), ansatz('b11'), ansatz('b02'), ansatz('c')
    lin_eqs = []
    for key, e in eqs:
        ex = e.subs(f, fexpr).subs({b20: A20, b11: A11, b02: A02, c: Ac}).doit()
        P = sp.Poly(canon(ex), S, C)
        lin_eqs += list(P.coeffs())
    unknowns = coeffs + list(params)
    M, rhs = sp.linear_eq_to_matrix(lin_eqs, unknowns)
    Mf = sp.lambdify((t2, t4), M, 'numpy')
    return Mf, len(coeffs), nparams

def residual(Mf, nb, nk, tv2, tv4, exclude=None):
    M = np.array(Mf(tv2, tv4), dtype=float)
    Mb, Mk = M[:, :nb], M[:, nb:]
    if exclude is not None:      # restrict Killing parameters to the complement of a trivial direction (H^2 for w=2)
        e = exclude/np.linalg.norm(exclude)
        Bc = np.linalg.svd(np.eye(len(e)) - np.outer(e, e))[0][:, :len(e)-1]
        Mk = Mk @ Bc
    Q, _ = np.linalg.qr(Mb)
    r = np.linalg.matrix_rank(Mb, tol=1e-9)
    Q = Q[:, :r]
    N = Mk - Q @ (Q.T @ Mk)
    sv = np.linalg.svd(N, compute_uv=False)
    return sv.min() if sv.size else 0.0

if __name__ == '__main__':
    w = int(sys.argv[1]); muval = sp.Rational(sys.argv[2]) if len(sys.argv) > 2 else 1
    Mf, nb, nk = build(w, muval=muval)
    excl = np.array([1.0, 2.0, 1.0]) if w == 2 else None
    grid = np.linspace(-1.0, 1.0, 81)
    res = np.zeros((len(grid), len(grid)))
    for i, a in enumerate(grid):
        for j, b in enumerate(grid):
            res[i, j] = residual(Mf, nb, nk, a, b, excl)
    np.save(f'scan_w{w}_mu{str(muval).replace("/","o")}.npy', res)
    thr = 1e-8
    hits = [(grid[i], grid[j], res[i, j]) for i in range(len(grid)) for j in range(len(grid)) if res[i, j] < thr]
    print(f"w={w}, mu={muval}: grid {len(grid)}x{len(grid)}, min residual {res.min():.2e}, hits below {thr}: {len(hits)}")
    for h in hits[:40]: print("   t2=%.4f t4=%.4f res=%.2e" % h)
    # also report near-minima (local minima below 1e-3) to catch curves
    from scipy.ndimage import minimum_filter
    loc = (res == minimum_filter(res, size=5)) & (res < 1e-3)
    print("  local minima below 1e-3:", [(round(grid[i],3), round(grid[j],3), '%.1e' % res[i,j]) for i, j in zip(*np.where(loc))][:30])
