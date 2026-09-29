import sympy as sp, sys, time, itertools, pickle
th, r = sp.symbols('theta r', positive=True)
px, pz = sp.symbols('p_x p_z')
S, C = sp.symbols('S C')
J = 10
Fj = sp.symbols(f'F0:{J}')

def dq(expr, which):
    if which == 'x': return sp.sin(th)*sp.diff(expr, r) + sp.cos(th)/r*sp.diff(expr, th)
    return sp.cos(th)*sp.diff(expr, r) - sp.sin(th)/r*sp.diff(expr, th)

def redSC(x):
    x = sp.expand(x); P = sp.Poly(x, C); out = 0
    for (k,), co in zip(P.monoms(), P.coeffs()): out += co*(1-S**2)**(k//2)*C**(k % 2)
    return sp.expand(out)

def canon(e):
    e = sp.expand(sp.expand_trig(e).subs({sp.tan(th): sp.sin(th)/sp.cos(th)}))
    e = sp.together(e); num, den = sp.fraction(e)
    rd = lambda x: redSC(sp.expand(x.subs({sp.sin(th): S, sp.cos(th): C})))
    return rd(num), rd(den)

def build(n, d, params, muval):
    """Return equations (jets) and unknown-function jet symbols for degree-n integral with A = L^d Q_{n-d}."""
    x = r*sp.sin(th); z = r*sp.cos(th); L = x*pz - z*px
    f = sp.Function('f')(th); V = r**muval*f
    # Killing tensor, even-even monomials only
    q = 0; k = 0
    for i in range(n-d, -1, -1):
        j = n-d-i
        if i % 2 == 0 and j % 2 == 0:
            q += params[k]*px**i*pz**j; k += 1
    assert k == len(params), (k, len(params))
    A = L**d*q
    unknowns = {}   # (level m, a) -> Function
    I = A
    for m in range(n-2, -1, -2):
        e = d + muval*(n-m)//2 if (n-m) % 2 == 0 else None
        e = d + muval*sp.Rational(n-m, 2)
        for a in range(m, -1, -1):
            u = sp.Function(f'u{m}_{a}')(th); unknowns[(m, a)] = u
            I += r**e*u*px**a*pz**(m-a)
    Vx, Vz = dq(V, 'x'), dq(V, 'z')
    PB = px*dq(I, 'x') + pz*dq(I, 'z') - (sp.diff(I, px)*Vx + sp.diff(I, pz)*Vz)
    PB = sp.expand(PB)
    P = sp.Poly(PB, px, pz)
    eqs = {}
    for (i, j), co in zip(P.monoms(), P.coeffs()):
        m = i + j
        if m == n + 1: 
            assert sp.simplify(co) == 0, "Killing tensor failure"
            continue
        co = sp.simplify(co)
        # remove the r power
        co = sp.powsimp(sp.expand(co), force=True)
        rp = [sp.powsimp(t).as_powers_dict().get(r, 0) for t in sp.Add.make_args(co)]
        assert len(set(rp)) == 1, (m, i, j, set(rp))
        co = sp.expand(sp.powsimp(co/r**rp[0], force=True))
        assert not co.has(r)
        eqs[(m, i, j)] = co
    # jets
    Uj = {key: sp.symbols(f'U{key[0]}_{key[1]}_0:{J}') for key in unknowns}
    rep = {}
    for key, u in unknowns.items():
        for k in range(J-1, -1, -1): rep[sp.Derivative(u, (th, k)) if k else u] = Uj[key][k]
    for k in range(J-1, -1, -1): rep[sp.Derivative(f, (th, k)) if k else f] = Fj[k]
    eqsJ = {key: sp.expand(e.doit().subs(rep)) for key, e in eqs.items()}
    return eqsJ, Uj, unknowns

def eliminate(n, d, params, muval, verbose=True):
    t0 = time.time()
    eqsJ, Uj, unknowns = build(n, d, params, muval)
    if verbose: print(f"  system built ({len(eqsJ)} eqs, {len(Uj)} unknown functions) in {time.time()-t0:.0f}s", flush=True)
    def Dth(e):
        out = sp.diff(e, th)
        for key in Uj:
            for k in range(J-1): out += sp.diff(e, Uj[key][k])*Uj[key][k+1]
        for k in range(J-1): out += sp.diff(e, Fj[k])*Fj[k+1]
        return sp.expand(out)
    # level by level: solve for first derivatives sequentially, collect algebraic relations
    subsP = {}; relations = []
    levels = sorted({key[0] for key in eqsJ}, reverse=True)
    for m in levels:
        if m == 1:
            # unknown D = u0_0: two equations cos D' + e sin D = T1 ; -sin D' + e cos D = T0
            D0, D1 = Uj[(0, 0)][0], Uj[(0, 0)][1]
            E1 = sp.expand(eqsJ[(1, 1, 0)].subs(subsP)); E0 = sp.expand(eqsJ[(1, 0, 1)].subs(subsP))
            wf = sp.simplify(E1.coeff(D0)/sp.sin(th))
            T1 = -sp.expand(E1.subs({D0: 0, D1: 0})); T0 = -sp.expand(E0.subs({D0: 0, D1: 0}))
            Dsol = sp.expand((sp.sin(th)*T1 + sp.cos(th)*T0)/wf); Dp = sp.expand(sp.cos(th)*T1 - sp.sin(th)*T0)
            relations.append(sp.expand(sp.numer(sp.together((Dth(Dsol) - Dp).subs(subsP)))))
            continue
        for a in range(m, 0, -1):
            e = sp.expand(eqsJ[(m, a, m-a)].subs(subsP))
            var = Uj[(m-1, a-1)][1]
            sol = sp.solve(e, var)
            assert len(sol) == 1, (m, a)
            subsP[var] = sp.expand(sol[0])
            # keep subsP closed: substitute into previously stored expressions
            for k2 in list(subsP): subsP[k2] = sp.expand(subsP[k2].subs(var, subsP[var]))
        e = sp.expand(eqsJ[(m, 0, m)].subs(subsP))
        relations.append(sp.numer(sp.together(e)))
    if verbose: print(f"  {len(relations)} algebraic relations, in {time.time()-t0:.0f}s", flush=True)
    u0 = [Uj[key][0] for key in Uj if key[0] > 0]      # unknown functions to eliminate (D already algebraic)
    def red(e): return sp.expand(e.subs(subsP))
    # generate derivatives of relations until we have >= len(u0)+ extra
    pool = []
    nu = len(u0); nrel = len(relations)
    # derivative counts: enough rows to solve for all unknowns plus a few extra conditions
    counts = [nu - (nrel - 1)] + [1]*(nrel - 1)      # rows = nu + ... ; then add one extra derivative each
    counts = [c + 1 for c in counts]
    for R, cnt in zip(relations, counts):
        cur = sp.expand(red(R))
        for k in range(cnt):
            pool.append(cur); cur = sp.expand(red(Dth(cur)))
    # canonical forms
    poolC = []
    for e in pool:
        nn, dd = canon(e); poolC.append(nn)
    M, rhs = sp.linear_eq_to_matrix(poolC, u0)
    if verbose: print(f"  jet matrix {M.shape}, built in {time.time()-t0:.0f}s", flush=True)
    # choose nu rows with nonzero determinant to solve, remaining rows give conditions
    rows = list(range(M.shape[0])); best = None
    for comb in itertools.combinations(rows, nu):
        sub = M.extract(list(comb), list(range(nu)))
        dt = sp.cancel(sub.det())
        if dt != 0:
            best = comb; detv = dt; break
    if best is None: raise RuntimeError("all subsystems singular")
    if verbose: print("  solving subsystem rows", best, "det:", sp.factor(redSC(sp.numer(sp.together(detv)))), flush=True)
    sub = M.extract(list(best), list(range(nu))); rsub = rhs.extract(list(best), [0])
    usol = sub.LUsolve(rsub); usol = [sp.cancel(sp.together(x)) for x in usol]
    conds = []
    for i in rows:
        if i in best: continue
        v = sp.cancel(sp.together((M[i, :]*sp.Matrix(usol))[0] - rhs[i]))
        v = sp.numer(v)
        conds.append(sp.factor(redSC(v)))
    if verbose:
        print(f"  {len(conds)} conditions in {time.time()-t0:.0f}s", flush=True)
        for c in conds: print("   COND:", str(c)[:900], flush=True)
    return M, rhs, u0, poolC, usol, conds

if __name__ == '__main__':
    n, d, muval = int(sys.argv[1]), int(sys.argv[2]), sp.Rational(sys.argv[3])
    nparams = len([1 for i in range(n-d, -1, -1) if i % 2 == 0 and (n-d-i) % 2 == 0])
    params = sp.symbols(f'a0:{nparams}')
    tag = ''
    if len(sys.argv) > 4:
        k = int(sys.argv[4]); params = tuple(sp.Integer(1) if i == k else sp.Integer(0) for i in range(nparams)); tag = f'_k{k}'
    M, rhs, u0, pool, usol, conds = eliminate(n, d, params, muval)
    pickle.dump((M, rhs, u0, pool, usol, conds, params), open(f'deg{n}_d{d}_mu{str(muval).replace("/","o")}{tag}.pkl', 'wb'))
    print("saved")
