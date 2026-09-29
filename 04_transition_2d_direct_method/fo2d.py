import sympy as sp, sys, time, itertools, pickle
r, th, pr, pt = sp.symbols('r theta p_r p_theta')
alpha, beta, delta = sp.symbols('alpha beta delta')
J = 8

def legendre_modes(l):
    """P_l(cos theta) as dict m -> coefficient of exp(i m theta)"""
    e = sp.expand(sp.legendre(l, sp.cos(th)).rewrite(sp.exp))
    e = sp.expand(e)
    modes = {}
    for m in range(-l, l+1):
        c = e.coeff(sp.exp(sp.I*m*th)) if m != 0 else e.subs(sp.exp(sp.I*th), 0).subs(sp.exp(-sp.I*th), 0)
        # robust extraction: substitute w = exp(i theta) and take polynomial coefficient
    w = sp.Symbol('w')
    ew = sp.expand(e.subs(sp.exp(sp.I*th), w).subs(sp.exp(-sp.I*th), 1/w))
    P = sp.Poly(sp.expand(ew*w**l), w)
    for (k,), c in zip(P.monoms(), P.coeffs()):
        modes[k - l] = sp.nsimplify(c)
    return modes

def mode_system(m, ls, A, Bjets):
    """Equations of the O(eps) condition in Fourier mode m: {Q_m e^{im th}, H0} = -{I0, H1}_m.
       Unknown Q_m = sum c_{jk}(r) p_r^j p_th^k (j+k even <= 4). Returns list of (equation) linear in c-jets and B-jets."""
    cs = {}
    Q = 0
    for j in range(5):
        for k in range(5 - j):
            if (j + k) % 2 == 0:
                cs[(j, k)] = sp.Function(f'c{j}{k}')(r); Q += cs[(j, k)]*pr**j*pt**k
    Ap = sp.diff(A, r)
    # {F, H0} with F = Q e^{im th}:  p_r dF/dr + (p_th/r^2) dF/dth - (A' - p_th^2/r^3) dF/dp_r
    lhs = pr*sp.diff(Q, r) + (pt/r**2)*(sp.I*m)*Q - (Ap - pt**2/r**3)*sp.diff(Q, pr)
    # source: -{I0, H1}_m where I0 = (alpha + beta H0) L^2 + delta L^4, H1 = sum_l B_l(r) P_l
    H0 = pr**2/2 + pt**2/(2*r**2) + A
    src = 0
    for l in ls:
        modes = legendre_modes(l)
        if m not in modes: continue
        cm = modes[m]
        Bl = Bjets[l][0]; Blp = Bjets[l][1]
        H1m = Bl*cm            # times e^{im th}
        dH1_dr = Blp*cm; dH1_dth = sp.I*m*Bl*cm
        # {I0, H1} = F'(H0) L^2 {H0,H1} + (F(H0) + 2 delta L^2) {L^2, H1};  {H0,H1} = p_r dH1/dr + (p_th/r^2) dH1/dth ; {L^2,H1} = -2 L dH1/dth
        F = alpha + beta*H0; Fp = beta
        br = Fp*pt**2*(pr*dH1_dr + (pt/r**2)*dH1_dth) + (F + 2*delta*pt**2)*(-2*pt*dH1_dth)
        src += -br
    eq = sp.expand(lhs - src)
    P = sp.Poly(eq, pr, pt)
    eqs = [(mon, sp.together(co)) for mon, co in zip(P.monoms(), P.coeffs())]
    return eqs, cs

def eliminate_mode(m, ls, A, D=3, verbose=True):
    t0 = time.time()
    Bjets = {l: sp.symbols(f'B{l}_0:{J}') for l in ls}
    eqs, cs = mode_system(m, ls, A, Bjets)
    keys = sorted(cs)
    Cj = {key: sp.symbols(f'C{key[0]}{key[1]}_0:{J}') for key in keys}
    rep = {}
    for key in keys:
        for k in range(J-1, -1, -1): rep[sp.Derivative(cs[key], (r, k)) if k else cs[key]] = Cj[key][k]
    def D_(e):
        out = sp.diff(e, r)
        for key in keys:
            for k in range(J-1): out += sp.diff(e, Cj[key][k])*Cj[key][k+1]
        for l in ls:
            for k in range(J-1): out += sp.diff(e, Bjets[l][k])*Bjets[l][k+1]
        return sp.expand(out)
    rows = []
    for mon, e in eqs:
        ej = sp.expand(sp.numer(sp.together(e.subs(rep).doit())))
        cur = ej
        for d in range(D+1):
            rows.append(cur); cur = sp.expand(sp.numer(sp.together(D_(cur))))
    cvars = [v for key in keys for v in Cj[key] if any(rw.has(v) for rw in rows)]
    M, rhs = sp.linear_eq_to_matrix(rows, cvars)
    if verbose: print(f"  mode {m}: {len(rows)} rows, {len(cvars)} c-jets, built in {time.time()-t0:.0f}s", flush=True)
    from sympy.polys.matrices import DomainMatrix
    K_ = sp.QQ_I.frac_field(r)
    Md = DomainMatrix.from_Matrix(M).convert_to(K_)
    LN = Md.transpose().nullspace().to_Matrix()      # rows are left null vectors
    conds = []
    for i in range(LN.shape[0]):
        n = LN[i, :]
        c = sp.expand(sp.numer(sp.together(sum(n[j]*rhs[j] for j in range(len(rhs))))))
        if c != 0: conds.append(sp.factor(c))
    if verbose:
        print(f"  mode {m}: {LN.shape[0]} left-null vectors, {len(conds)} nontrivial conditions in {time.time()-t0:.0f}s", flush=True)
        for c in conds[:8]: print("   COND:", str(c)[:700], flush=True)
    return conds, Bjets

if __name__ == '__main__':
    A = -1/(1+r)
    m = int(sys.argv[1]); ls = [int(x) for x in sys.argv[2].split(',')]
    conds, Bjets = eliminate_mode(m, ls, A)
    pickle.dump((conds, [str(b) for l in ls for b in Bjets[l]]), open(f'fo2d_m{m}.pkl', 'wb'))
