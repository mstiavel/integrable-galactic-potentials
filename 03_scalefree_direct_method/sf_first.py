import sympy as sp, sys, time, pickle, random
from sf_general import build, Fj, J, S, C, redSC, th
from sf_fast import to_SC, Dsc

def first_condition(n, d, muval, kunit=None, nder=None, verbose=True):
    t0 = time.time()
    nparams = len([1 for i in range(n-d, -1, -1) if i % 2 == 0 and (n-d-i) % 2 == 0])
    params = sp.symbols(f'a0:{nparams}')
    if kunit is not None: params = tuple(sp.Integer(1) if i == kunit else sp.Integer(0) for i in range(nparams))
    eqsJ, Uj, unknowns = build(n, d, params, muval)
    eqs = {}
    for key, e in eqsJ.items():
        num, den = to_SC(e); eqs[key] = sp.expand(num)
    m = n - 1
    def D(e):
        out = Dsc(e)
        for key in Uj:
            for k in range(J-1): out += sp.diff(e, Uj[key][k])*Uj[key][k+1]
        for k in range(J-1): out += sp.diff(e, Fj[k])*Fj[k+1]
        return sp.expand(out)
    subsP = {}
    def apply(e):
        e = sp.expand(e.subs(subsP)); e = sp.together(e); num, den = sp.fraction(e)
        return redSC(sp.expand(num))
    for a in range(m, 0, -1):
        e = apply(eqs[(m, a, m-a)])
        var = Uj[(m-1, a-1)][1]
        co = sp.Poly(e, var).coeff_monomial(var); rest = sp.expand(e - co*var)
        subsP[var] = sp.cancel(-rest/co)
        for k2 in list(subsP): subsP[k2] = sp.cancel(subsP[k2].subs(var, subsP[var]))
    R1 = apply(eqs[(m, 0, m)])
    u0 = [Uj[(m-1, a)][0] for a in range(m)]
    nu = len(u0)
    if nder is None: nder = nu
    pool = []; cur = R1
    for k in range(nder + 1):
        pool.append(cur); cur = apply(D(cur))
    M, rhs = sp.linear_eq_to_matrix(pool, u0)
    if verbose: print(f"  level-{m} relation and {nder} derivatives: matrix {M.shape} in {time.time()-t0:.0f}s", flush=True)
    rng = random.Random(5)
    num = {s_: sp.Rational(rng.randint(2, 40), rng.randint(2, 40)) for s_ in list(Fj) + [S, C]}
    Mn = M.subs(num)
    cols = []
    for j in range(nu):
        if Mn.extract(list(range(M.shape[0])), cols + [j]).rank() > len(cols): cols.append(j)
    rows = []
    for i in range(M.shape[0]):
        if Mn.extract(rows + [i], cols).rank() > len(rows): rows.append(i)
    free = [u0[j] for j in range(nu) if j not in cols]
    if verbose: print("  rank", len(rows), "free:", free, flush=True)
    from sympy.polys.matrices import DomainMatrix
    sub = M.extract(rows, cols)
    rsub = rhs.extract(rows, [0]) - (M.extract(rows, [j for j in range(nu) if j not in cols])*sp.Matrix(free) if free else sp.zeros(len(rows), 1))
    rest = [i for i in range(M.shape[0]) if i not in rows]
    allsyms = set()
    for e in list(sub) + list(rsub) + list(M.extract(rest, list(range(nu)))) + list(rhs): allsyms |= e.free_symbols
    K_ = sp.QQ.frac_field(*sorted(allsyms, key=str))
    A = DomainMatrix.from_Matrix(sub).convert_to(K_); B = DomainMatrix.from_Matrix(rsub).convert_to(K_)
    X = A.lu_solve(B)
    Mrest = DomainMatrix.from_Matrix(M.extract(rest, cols)).convert_to(K_)
    rrest = rhs.extract(rest, [0]) - (M.extract(rest, [j for j in range(nu) if j not in cols])*sp.Matrix(free) if free else sp.zeros(len(rest), 1))
    Res = (Mrest*X - DomainMatrix.from_Matrix(rrest).convert_to(K_)).to_Matrix()
    conds = [sp.factor(redSC(sp.expand(sp.numer(v)))) for v in Res]
    if verbose:
        print(f"  {len(conds)} conditions in {time.time()-t0:.0f}s", flush=True)
        for c in conds: print("   COND:", str(c)[:1500], flush=True)
    return conds, params

if __name__ == '__main__':
    n, d, muval = int(sys.argv[1]), int(sys.argv[2]), sp.Rational(sys.argv[3])
    kunit = int(sys.argv[4]) if len(sys.argv) > 4 else None
    conds, params = first_condition(n, d, muval, kunit)
    pickle.dump((conds, params), open(f'first_deg{n}_d{d}' + (f'_k{kunit}' if kunit is not None else '') + '.pkl', 'wb'))
