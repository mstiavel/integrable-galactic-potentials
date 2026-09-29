import sympy as sp, sys, time, pickle
from fo2d import legendre_modes
r, th, pr, pt = sp.symbols('r theta p_r p_theta')
Lg = sp.Symbol('Lg')                      # stands for log(r/(1+r)) when the seed is Jaffe
J = 9

SEEDS = {
    'hernquist': (-1/(1+r), 1/(1+r)**2, None),
    'jaffe':     (Lg, 1/(r*(1+r)), 1/(r*(1+r))),          # A = log(r/(1+r)), A' given, dLg/dr given
    'cored':     (-1/(1+r**2), 2*r/(1+r**2)**2, None),
}

def seed_integral(n, H0, L2):
    if n == 4:
        a, b, d = sp.symbols('alpha beta delta')
        return (a + b*H0)*L2 + d*L2**2, (a, b, d)
    if n == 6:
        a, b, g, d, e, k = sp.symbols('alpha beta gamma delta eta kappa')
        return (a + b*H0 + g*H0**2)*L2 + (d + e*H0)*L2**2 + k*L2**3, (a, b, g, d, e, k)

def mode_system(m, ls, n, A, Ap, Bjets):
    cs = {}; Q = 0
    for j in range(n+1):
        for k in range(n+1-j):
            if (j + k) % 2 == 0:
                cs[(j, k)] = sp.Function(f'c{j}{k}')(r); Q += cs[(j, k)]*pr**j*pt**k
    lhs = pr*sp.diff(Q, r) + (pt/r**2)*(sp.I*m)*Q - (Ap - pt**2/r**3)*sp.diff(Q, pr)
    H0 = pr**2/2 + pt**2/(2*r**2) + A
    I0, params = seed_integral(n, H0, pt**2)
    # {I0, H1} = dI0/dpr * (-dH1/dr) ... use full bracket: {I0,H1} = I0_r H1_pr - I0_pr H1_r + I0_th H1_pt - I0_pt H1_th, H1 = B(r) P e^{im th}
    src = 0
    for l in ls:
        modes = legendre_modes(l)
        if m not in modes: continue
        cm = modes[m]; Bl, Blp = Bjets[l][0], Bjets[l][1]
        H1r = Blp*cm; H1th = sp.I*m*Bl*cm
        br = -sp.diff(I0, pr)*H1r - sp.diff(I0, pt)*H1th
        src += -br
    eq = sp.expand(lhs - src)
    P = sp.Poly(eq, pr, pt)
    return [(mon, sp.together(co)) for mon, co in zip(P.monoms(), P.coeffs())], cs, params

def eliminate_mode(m, ls, n, seed, D=None, verbose=True):
    t0 = time.time()
    A, Ap, dLg = SEEDS[seed]
    Bjets = {l: sp.symbols(f'B{l}_0:{J}') for l in ls}
    eqs, cs, params = mode_system(m, ls, n, A, Ap, Bjets)
    keys = sorted(cs); Cj = {key: sp.symbols(f'C{key[0]}{key[1]}_0:{J}') for key in keys}
    rep = {}
    for key in keys:
        for k in range(J-1, -1, -1): rep[sp.Derivative(cs[key], (r, k)) if k else cs[key]] = Cj[key][k]
    def D_(e):
        out = sp.diff(e, r)
        if dLg is not None: out += sp.diff(e, Lg)*dLg
        for key in keys:
            for k in range(J-1): out += sp.diff(e, Cj[key][k])*Cj[key][k+1]
        for l in ls:
            for k in range(J-1): out += sp.diff(e, Bjets[l][k])*Bjets[l][k+1]
        return sp.expand(out)
    nunk = len(keys); neq = len(eqs)
    if D is None:
        D = 0
        while neq*(D+1) <= nunk*(D+2) + 2: D += 1
    rows = []
    for mon, e in eqs:
        cur = sp.expand(sp.numer(sp.together(e.subs(rep).doit())))
        for d in range(D+1):
            rows.append(cur); cur = sp.expand(sp.numer(sp.together(D_(cur))))
    cvars = [v for key in keys for v in Cj[key] if any(rw.has(v) for rw in rows)]
    M, rhs = sp.linear_eq_to_matrix(rows, cvars)
    if verbose: print(f"  seed={seed} n={n} mode {m}: {len(rows)} rows, {len(cvars)} c-jets (D={D}), built in {time.time()-t0:.0f}s", flush=True)
    from sympy.polys.matrices import DomainMatrix
    import random
    gens = [r] + ([Lg] if dLg is not None else [])
    K_ = sp.QQ_I.frac_field(*gens)
    # rank-aware square subsystem chosen numerically, then exact solve; conditions from the remaining rows
    rng = random.Random(7)
    num = {r: sp.Rational(rng.randint(2, 40), rng.randint(2, 40))}
    if dLg is not None: num[Lg] = sp.Rational(rng.randint(2, 40), rng.randint(2, 40))
    Mn = M.subs(num)
    cols = []
    for j in range(M.shape[1]):
        if Mn.extract(list(range(M.shape[0])), cols + [j]).rank() > len(cols): cols.append(j)
    rows_ = []
    for i in range(M.shape[0]):
        if Mn.extract(rows_ + [i], cols).rank() > len(rows_): rows_.append(i)
    free = [cvars[j] for j in range(M.shape[1]) if j not in cols]
    if verbose: print(f"  rank {len(rows_)} of {M.shape[1]} jets; {len(free)} free jets", flush=True)
    sub = M.extract(rows_, cols)
    rsub = rhs.extract(rows_, [0]) - (M.extract(rows_, [j for j in range(M.shape[1]) if j not in cols])*sp.Matrix(free) if free else sp.zeros(len(rows_), 1))
    # rhs is linear in (param x B-jet) monomials and in the free jets: decompose and solve multi-column over K_
    rest = [i for i in range(M.shape[0]) if i not in rows_]
    rrest = rhs.extract(rest, [0]) - (M.extract(rest, [j for j in range(M.shape[1]) if j not in cols])*sp.Matrix(free) if free else sp.zeros(len(rest), 1))
    lin_syms = list(params) + [sy for l in ls for sy in Bjets[l]] + list(free)
    monos = set()
    for e in list(rsub) + list(rrest):
        P_ = sp.Poly(sp.expand(e), *lin_syms); monos |= set(P_.monoms())
    monos = sorted(monos)
    def col(vec, mono):
        return sp.Matrix([sp.Poly(sp.expand(e), *lin_syms).coeff_monomial(sp.prod([sy**k for sy, k in zip(lin_syms, mono)])) for e in vec])
    Bmat = sp.Matrix.hstack(*[col(rsub, mo) for mo in monos])
    Rmat = sp.Matrix.hstack(*[col(rrest, mo) for mo in monos])
    Xall = DomainMatrix.from_Matrix(sub).convert_to(K_).lu_solve(DomainMatrix.from_Matrix(Bmat).convert_to(K_))
    Mrest = DomainMatrix.from_Matrix(M.extract(rest, cols)).convert_to(K_)
    Res = (Mrest*Xall - DomainMatrix.from_Matrix(Rmat).convert_to(K_)).to_Matrix()
    Res = sp.Matrix([[sum(Res[i, k]*sp.prod([sy**e for sy, e in zip(lin_syms, monos[k])]) for k in range(len(monos)))] for i in range(len(rest))])
    class _LN: shape = (len(rest), 0)
    LN = _LN()
    conds = []
    for v in Res:
        c = sp.expand(sp.numer(sp.together(v)))
        if c != 0:
            if any(c.has(fu) for fu in free): print("   WARNING: condition involves a free jet", flush=True)
            conds.append(c)
    if verbose: print(f"  {len(conds)} nontrivial conditions in {time.time()-t0:.0f}s", flush=True)
    return conds, Bjets, params

def analyse(conds, Bjets, l, params, dLg=None):
    """Joint solution space for a single multipole l: jet nullspace and consistency; report B'/B or the ODE."""
    B = Bjets[l]
    core = []
    for c in conds:
        fl = sp.factor_list(c)[1]
        fs = [f for f, mlt in fl if f.has(B[0])]
        if fs: core.append(sp.expand(fs[0]))
    if not core: print("   no conditions on B"); return
    kmax = max(max(k for k in range(J) if c.has(B[k])) for c in core)
    M, rhs = sp.linear_eq_to_matrix(core, list(B[:kmax+1]))
    from sympy.polys.matrices import DomainMatrix
    gens = [r] + ([Lg] if dLg is not None else []) + list(params)
    K_ = sp.QQ.frac_field(*gens)
    ns = DomainMatrix.from_Matrix(M).convert_to(K_).nullspace().to_Matrix()
    dim = ns.shape[0]
    print(f"   {len(core)} conditions (max order {kmax}); jet-nullspace dimension {dim}")
    def Dr(e):
        return sp.diff(e, r) + (sp.diff(e, Lg)*dLg if dLg is not None else 0)
    if dim == 1:
        v = [sp.cancel(x) for x in ns[0, :]]
        phi = [sp.cancel(v[k]/v[0]) for k in range(kmax+1)]
        print("   B'/B =", sp.factor(phi[1]))
        cur = phi[1]
        for k in range(2, kmax+1):
            nxt = sp.cancel(Dr(cur) + phi[1]*cur)
            c = sp.factor(sp.numer(sp.cancel(nxt - phi[k])))
            print(f"   consistency at order {k}:", 0 if c == 0 else str(c)[:200])
            cur = phi[k]
    elif dim == 2:
        # second-order ODE: express B'' via B, B'
        sub = M[:, 2:kmax+1]; rr = -(M[:, 0:2]*sp.Matrix(B[:2]))
        Kb = sp.QQ.frac_field(*(gens + list(B[:2])))
        X = DomainMatrix.from_Matrix(sub).convert_to(Kb).lu_solve(DomainMatrix.from_Matrix(rr).convert_to(Kb)).to_Matrix()
        B2e = sp.cancel(X[0])
        a = sp.cancel(B2e.subs(B[1], 0)/B[0]); b = sp.cancel(B2e.subs(B[0], 0)/B[1])
        print("   B'' = a B + b B'  with a =", sp.factor(a), "; b =", sp.factor(b))
        s = sp.Symbol('s')
        ind = sp.factor(sp.limit(sp.expand(s*(s-1) - b*r*s - a*r**2), r, 0))
        print("   indicial polynomial at r=0:", ind)
    else:
        print("   solution space of dimension", dim)

if __name__ == '__main__':
    seed, n, m = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); ls = [int(x) for x in sys.argv[4].split(',')]
    conds, Bjets, params = eliminate_mode(m, ls, n, seed)
    pickle.dump((conds, params), open(f'fo2d_{seed}_n{n}_m{m}.pkl', 'wb'))
    analyse(conds, Bjets, max(ls), params, SEEDS[seed][2])
