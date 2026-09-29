import sympy as sp, pickle, sys, time
r, th, pr, pt = sp.symbols('r theta p_r p_theta')
alpha, beta, delta, b = sp.symbols('alpha beta delta b')
A = -1/(1+r**2); Ap = sp.diff(A, r)
B2 = r**2/(1+r**2)**2
H0 = pr**2/2 + pt**2/(2*r**2) + A
L2 = pt**2
I0 = (alpha + beta*H0)*L2 + delta*L2**2
# P2 = 1/4 + (3/8)(e^{2i th} + e^{-2i th})
modes_P2 = {0: sp.Rational(1, 4), 2: sp.Rational(3, 8), -2: sp.Rational(3, 8)}

def bracket_I0_H1_mode(m):
    """{I0, H1}_m for H1 = b B2 P2 (coefficient of e^{im th})"""
    cm = modes_P2[m]
    H1r = b*sp.diff(B2, r)*cm; H1th = sp.I*m*b*B2*cm
    return sp.expand(-sp.diff(I0, pr)*H1r - sp.diff(I0, pt)*H1th)

def solve_mode(m, degnum=14, powr=6, powq=7):
    """solve {Q_m e^{im th}, H0} = -{I0,H1}_m for Q_m quartic in momenta with rational coefficients"""
    t0 = time.time()
    unknown = []; Q = 0
    for j in range(5):
        for k in range(5 - j):
            if (j + k) % 2 == 0:
                cf = 0
                for ti, tf in enumerate([1, sp.atan(r), sp.log(1+r**2), sp.log(r), sp.atan(r)**2, sp.atan(r)*sp.log(1+r**2), sp.log(1+r**2)**2, sp.log(r)**2, sp.atan(r)*sp.log(r), sp.log(r)*sp.log(1+r**2)]):
                    for a in range(degnum + 1):
                        s = sp.Symbol(f'q{m}_{j}{k}_{ti}_{a}'); unknown.append(s); cf += s*r**a*tf
                Q += cf/(r**powr*(1+r**2)**powq)*pr**j*pt**k
    lhs = pr*sp.diff(Q, r) + (pt/r**2)*(sp.I*m)*Q - (Ap - pt**2/r**3)*sp.diff(Q, pr)
    T1, T2, T3 = sp.symbols('T1 T2 T3')
    eq = sp.expand(sp.numer(sp.together(lhs + bracket_I0_H1_mode(m))))
    eq = sp.expand(eq.subs({sp.atan(r): T1, sp.log(1+r**2): T2, sp.log(r): T3}))
    P = sp.Poly(eq, pr, pt, r, T1, T2, T3)
    lin = list(P.coeffs())
    sol = sp.solve(lin, unknown, dict=True)
    if not sol: print(f"  mode {m}: NO rational solution with this ansatz", flush=True); return None, unknown
    sol = sol[0]
    Qs = sp.simplify(Q.subs(sol))
    free = sorted(Qs.free_symbols - {r, pr, pt, alpha, beta, delta, b}, key=str)
    # check
    resid = sp.simplify(pr*sp.diff(Qs, r) + (pt/r**2)*(sp.I*m)*Qs - (Ap - pt**2/r**3)*sp.diff(Qs, pr) + bracket_I0_H1_mode(m))
    print(f"  mode {m}: solved in {time.time()-t0:.0f}s; residual {resid}; free constants {free}", flush=True)
    return Qs, free

if __name__ == '__main__':
    Q = {}
    for m in (2, 0):
        Qs, free = solve_mode(m)
        Q[m] = Qs
    # mode -2 is the conjugate-symmetric partner (real integral): Q_{-2} = Q_2 with i -> -i
    Q[-2] = Q[2].subs(sp.I, -sp.I) if Q[2] is not None else None
    pickle.dump(Q, open('cored_I1.pkl', 'wb'))
    for m in (2, 0):
        if Q[m] is not None:
            print(f"Q_{m} =", sp.factor(Q[m]))
