import sympy as sp, sys, time, pickle
from sf_direct import system, th, f
J = 8
Fj = sp.symbols(f'F0:{J}'); Bj = {n: sp.symbols(f'{n}_0:{J}') for n in ('b20','b11','b02','c')}
S, C = sp.symbols('S C')
muval = sp.Rational(sys.argv[1]); w = 4; params = sp.symbols('a b'); a, b = params
t0 = time.time()
eqs, fns = system(w, params, muval)
def jets(eqs):
    rep = {}
    for fn, name in zip(fns, ('b20','b11','b02','c')):
        for k in range(J-1, -1, -1): rep[sp.Derivative(fn, (th, k)) if k else fn] = Bj[name][k]
    for k in range(J-1, -1, -1): rep[sp.Derivative(f, (th, k)) if k else f] = Fj[k]
    return {key: sp.expand(e.doit().subs(rep)) for key, e in eqs}
E = jets(eqs)
def Dth(e):
    out = sp.diff(e, th)
    for name in Bj:
        for k in range(J-1): out += sp.diff(e, Bj[name][k])*Bj[name][k+1]
    for k in range(J-1): out += sp.diff(e, Fj[k])*Fj[k+1]
    return sp.expand(out)
E1, E2, E3, E4, E5, E6 = E[(3,3,0)], E[(3,2,1)], E[(3,1,2)], E[(3,0,3)], E[(1,1,0)], E[(1,0,1)]
b20p = sp.solve(E1, Bj['b20'][1])[0]; b02p = sp.solve(E4, Bj['b02'][1])[0]
b11p = sp.solve(E2.subs(Bj['b20'][1], b20p), Bj['b11'][1])[0]
subsP = {Bj['b20'][1]: b20p, Bj['b02'][1]: b02p, Bj['b11'][1]: b11p}
red = lambda e: sp.expand(e.subs(subsP))
R1 = sp.numer(sp.together(red(E3)))
c0, c1 = Bj['c'][0], Bj['c'][1]
wf = sp.simplify(E5.coeff(c0)/sp.sin(th))      # weight factor of c
T5 = -sp.expand(E5.subs({c0: 0, c1: 0})); T6 = -sp.expand(E6.subs({c0: 0, c1: 0}))
assert sp.simplify(E5 - (sp.cos(th)*c1 + wf*sp.sin(th)*c0 - T5)) == 0 and sp.simplify(E6 - (-sp.sin(th)*c1 + wf*sp.cos(th)*c0 - T6)) == 0
csol = sp.expand((sp.sin(th)*T5 + sp.cos(th)*T6)/wf); cp = sp.expand(sp.cos(th)*T5 - sp.sin(th)*T6)
R3 = sp.expand(red(sp.numer(sp.together(red(Dth(csol) - cp)))))
R1p = sp.expand(red(Dth(R1)))
Bvars = [Bj['b20'][0], Bj['b11'][0], Bj['b02'][0]]
M, rhs = sp.linear_eq_to_matrix([R1, R3, R1p], Bvars)
def canon(e):
    e = sp.expand(sp.expand_trig(e).subs({sp.tan(th): sp.sin(th)/sp.cos(th)})); e = sp.together(e); num, den = sp.fraction(e)
    def rd(x):
        x = sp.expand(x.subs({sp.sin(th): S, sp.cos(th): C})); P = sp.Poly(x, C); out = 0
        for (k,), co in zip(P.monoms(), P.coeffs()): out += co*(1-S**2)**(k//2)*C**(k % 2)
        return sp.expand(out)
    return rd(num), rd(den)
def redSC(x):
    x = sp.expand(x); P = sp.Poly(x, C); out = 0
    for (k,), co in zip(P.monoms(), P.coeffs()): out += co*(1-S**2)**(k//2)*C**(k % 2)
    return sp.expand(out)
Mn = sp.zeros(3,3); rn = sp.zeros(3,1)
for i in range(3):
    for j in range(3): n, d = canon(M[i,j]); Mn[i,j] = n/d
    n, d = canon(rhs[i]); rn[i] = n/d
dn, dd = canon(M.det()); print(f"mu={muval}: det =", sp.factor(dn), flush=True)
Bsol = [sp.cancel(sp.together(x)) for x in Mn.LUsolve(rn)]; subsB = dict(zip(Bvars, Bsol))
def DSC(e):
    out = sp.diff(e, S)*C - sp.diff(e, C)*S
    for k in range(J-1): out += sp.diff(e, Fj[k])*Fj[k+1]
    return out
Q = {}
for name, bs in zip(('b20','b11','b02'), Bsol):
    n, d = canon(subsP[Bj[name][1]]); bp = sp.cancel((n/d).subs(subsB))
    Q[name] = sp.numer(sp.cancel(sp.together(DSC(bs) - bp)))
al = {n: redSC(sp.Poly(q, Fj[3]).coeff_monomial(Fj[3])) for n, q in Q.items()}
be = {n: redSC(sp.Poly(q, Fj[3]).coeff_monomial(1)) for n, q in Q.items()}
for n in Q: print(n, "deg F3:", sp.Poly(Q[n], Fj[3]).degree(), " alpha:", sp.factor(al[n]), flush=True)
names = list(Q)
rel = sp.factor(redSC(al[names[0]]*be[names[1]] - al[names[1]]*be[names[0]]))
print("relation:", rel, flush=True)
fac = [fa for fa, m in sp.factor_list(rel)[1] if fa.has(Fj[1]) and (fa.has(a) or fa.has(b))]
master = fac[0]; print("MASTER ODE factor:", master, flush=True)
F2sol = sp.solve(master, Fj[2])[0]
D2 = lambda e: sp.diff(e, S)*C - sp.diff(e, C)*S + sp.diff(e, Fj[0])*Fj[1] + sp.diff(e, Fj[1])*Fj[2]
F3sol = sp.cancel(D2(F2sol).subs(Fj[2], F2sol))
for n, q in Q.items():
    v = sp.factor(redSC(sp.numer(sp.cancel(sp.together(q.subs(Fj[3], F3sol).subs(Fj[2], F2sol))))))
    print("residual", n, ":", v, flush=True)
pickle.dump((master, Q), open(f'mu_{str(muval).replace("/","o")}_w4.pkl','wb'))
print("time", time.time()-t0)
