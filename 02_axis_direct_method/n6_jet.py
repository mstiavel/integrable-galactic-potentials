import sympy as sp, itertools, sys, time
from n6_local2 import first_order_eqs
z, s = sp.symbols('z s')

def eliminate_jets(n, seedcoefs, A, D):
    eqs, gam, B1 = first_order_eqs(n, seedcoefs, A)
    eqs = [e for m, e in eqs if e != 0]
    gkeys = sorted(gam.keys())
    J = 12
    gj = {k: sp.symbols(f'G{k}_0:{J}') for k in gkeys}; bj = sp.symbols(f'Bj0:{J}')
    def tojet(e):
        rep = {}
        for k in gkeys:
            for j in range(J-1, -1, -1): rep[sp.Derivative(gam[k], (z, j)) if j else gam[k]] = gj[k][j]
        for j in range(J-1, -1, -1): rep[sp.Derivative(B1, (z, j)) if j else B1] = bj[j]
        return sp.expand(e.subs(rep))
    def Dz(e):
        out = sp.diff(e, z)
        for k in gkeys:
            for j in range(J-1): out += sp.diff(e, gj[k][j])*gj[k][j+1]
        for j in range(J-1): out += sp.diff(e, bj[j])*bj[j+1]
        return sp.expand(out)
    rows = []
    for e in eqs:
        je = tojet(e)
        cur = je
        for d in range(D+1):
            rows.append(cur); cur = Dz(cur)
    gvars = [v for k in gkeys for v in gj[k]]
    gvars = [v for v in gvars if any(r.has(v) for r in rows)]
    Mg = sp.Matrix([[r.coeff(v) for v in gvars] for r in rows])
    # left nullspace over Q(z): vectors c with c^T Mg = 0
    LN = Mg.T.nullspace()
    if not LN: return None
    bexprs = []
    for c in LN:
        comb = sum(ci*r for ci, r in zip(c, rows))
        comb = sp.expand(comb.subs({v: 0 for v in gvars}))
        comb = sp.numer(sp.together(comb))
        if comb != 0: bexprs.append(sp.expand(comb))
    return bexprs, bj

def indicial(expr, bj):
    order = max(j for j in range(12) if expr.has(bj[j]))
    terms = []
    for j in range(order+1):
        c = sp.expand(expr.coeff(bj[j]))
        if c == 0: continue
        P = sp.Poly(c, z); low = min(m[0] for m in P.monoms())
        terms.append((low - j, P.coeff_monomial(z**low)*sp.ff(s, j)))
    mn = min(t[0] for t in terms)
    return order, sp.factor(sum(c for m, c in terms if m == mn))

if __name__ == '__main__':
    a, b, d = sp.symbols('alpha beta delta')
    n = int(sys.argv[1]); D = int(sys.argv[2])
    seed = [a, b] if n == 4 else [a, b, d]
    t0 = time.time()
    res = eliminate_jets(n, seed, z, D)
    if res is None: print("no elimination at this D"); sys.exit()
    bexprs, bj = res
    print(f"n={n}, D={D}: {len(bexprs)} B1-only relations found in {time.time()-t0:.0f}s")
    for e in bexprs[:3]:
        order, ind = indicial(e, bj)
        print(f"   ODE order {order}, indicial at cusp: {ind}")
