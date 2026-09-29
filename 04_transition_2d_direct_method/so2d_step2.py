import sympy as sp, pickle, sys, time
from fo2d import legendre_modes
r, th, pr, pt = sp.symbols('r theta p_r p_theta')
alpha, beta, delta, b = sp.symbols('alpha beta delta b')
A = -1/(1+r**2); Ap = sp.diff(A, r); B2 = r**2/(1+r**2)**2
H0 = pr**2/2 + pt**2/(2*r**2) + A; I0 = (alpha + beta*H0)*pt**2 + delta*pt**4
J = 8
rp = sp.Symbol('r', positive=True)

def load_I1():
    s2 = pickle.load(open('cored_mode2_seq.pkl', 'rb')); s0 = pickle.load(open('cored_mode0_seq.pkl', 'rb'))
    fix2 = {sp.Symbol('K_c02'): 3*b*delta/2, sp.Symbol('K_c22'): 3*b*delta/4, sp.Symbol('K_c11'): 0, sp.Symbol('K_c00'): 3*alpha*b/4, sp.Symbol('K_c20'): 3*b*(alpha+beta)/8}
    fix0 = {sp.Symbol(f'K_c{n}'): 0 for n in ('40','31','22','13','04','20','11','02','00')}
    Q2 = 0; Q0 = 0
    for name, val in s2.items():
        j, k = int(name[1]), int(name[2]); Q2 += sp.simplify(val.subs(fix2)).subs(rp, r)*pr**j*pt**k
    for name, val in s0.items():
        j, k = int(name[1]), int(name[2]); Q0 += sp.simplify(val.subs(fix0)).subs(rp, r)*pr**j*pt**k
    Qm2 = Q2.subs(sp.I, -sp.I)
    return sp.simplify(Q2), sp.simplify(Q0), sp.simplify(Qm2)

def S2_modes(Q2, Q0, Qm2):
    """{I1, H1} with I1 = Q2 e^{2i th} + Q0 + Qm2 e^{-2i th}, H1 = b B2 (1/4 + 3/8 e^{2i th} + 3/8 e^{-2i th}); returns dict mode -> expr"""
    I1 = {2: Q2, 0: Q0, -2: Qm2}; H1 = {0: b*B2/4, 2: sp.Rational(3, 8)*b*B2, -2: sp.Rational(3, 8)*b*B2}
    out = {}
    for m1, Q in I1.items():
        for m2, h in H1.items():
            # {Q e^{i m1 th}, h e^{i m2 th}} = e^{i(m1+m2) th} [ Q_r h_pr - Q_pr h_r + (i m1 Q) h_pt - Q_pt (i m2 h) ]
            br = -sp.diff(Q, pr)*sp.diff(h, r) - sp.diff(Q, pt)*(sp.I*m2*h)
            out[m1+m2] = out.get(m1+m2, 0) + br
    return {m: sp.expand(e) for m, e in out.items()}

def eliminate_mode2(m, S2m, ls=(0, 2, 4), D=3, verbose=True):
    t0 = time.time()
    Cj = {l: sp.symbols(f'C{l}_0:{J}') for l in ls}
    cs = {}; Q = 0
    for j in range(5):
        for k in range(5 - j):
            if (j + k) % 2 == 0:
                cs[(j, k)] = sp.Function(f'c{j}{k}')(r); Q += cs[(j, k)]*pr**j*pt**k
    lhs = pr*sp.diff(Q, r) + (pt/r**2)*(sp.I*m)*Q - (Ap - pt**2/r**3)*sp.diff(Q, pr)
    # {I0, H2}_m
    src = 0
    for l in ls:
        modes = legendre_modes(l)
        if m not in modes: continue
        cm = modes[m]; Cl, Clp = Cj[l][0], Cj[l][1]
        src += -(-sp.diff(I0, pr)*Clp*cm - sp.diff(I0, pt)*(sp.I*m*Cl*cm))
    eq = sp.expand(lhs + S2m - src)            # {Q e^{im th},H0} + {I1,H1}_m + {I0,H2}_m = 0
    P = sp.Poly(eq, pr, pt)
    eqs = [sp.together(co) for co in P.coeffs()]
    keys = sorted(cs); cj = {key: sp.symbols(f'c{key[0]}{key[1]}_0:{J}') for key in keys}
    rep = {}
    for key in keys:
        for k in range(J-1, -1, -1): rep[sp.Derivative(cs[key], (r, k)) if k else cs[key]] = cj[key][k]
    def D_(e):
        out = sp.diff(e, r)
        for key in keys:
            for k in range(J-1): out += sp.diff(e, cj[key][k])*cj[key][k+1]
        for l in ls:
            for k in range(J-1): out += sp.diff(e, Cj[l][k])*Cj[l][k+1]
        return sp.expand(out)
    rows = []
    for e in eqs:
        cur = sp.expand(sp.numer(sp.together(e.subs(rep).doit())))
        for d in range(D+1):
            rows.append(cur); cur = sp.expand(sp.numer(sp.together(D_(cur))))
    cvars = [v for key in keys for v in cj[key] if any(rw.has(v) for rw in rows)]
    M, rhs = sp.linear_eq_to_matrix(rows, cvars)
    if verbose: print(f"  mode {m}: {len(rows)} rows, {len(cvars)} jets, built in {time.time()-t0:.0f}s", flush=True)
    from sympy.polys.matrices import DomainMatrix
    K_ = sp.QQ_I.frac_field(r)
    LN = DomainMatrix.from_Matrix(M).convert_to(K_).transpose().nullspace().to_Matrix()
    conds = []
    for i in range(LN.shape[0]):
        c = sp.expand(sp.numer(sp.together(sum(LN[i, j]*rhs[j] for j in range(len(rhs))))))
        if c != 0: conds.append(c)
    if verbose: print(f"  mode {m}: {LN.shape[0]} left-null vectors, {len(conds)} conditions in {time.time()-t0:.0f}s", flush=True)
    return conds, Cj

if __name__ == '__main__':
    m = int(sys.argv[1])
    Q2, Q0, Qm2 = load_I1()
    S2 = S2_modes(Q2, Q0, Qm2)
    print("S2 modes present:", {k: (0 if v == 0 else 'nonzero') for k, v in S2.items()}, flush=True)
    conds, Cj = eliminate_mode2(m, S2.get(m, 0))
    pickle.dump(conds, open(f'so2d_m{m}_conds.pkl', 'wb'))
    for c in conds: print("  COND:", str(sp.factor(c))[:600], flush=True)
