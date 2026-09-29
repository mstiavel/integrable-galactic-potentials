import sympy as sp, sys, time
from sf_classify import classify, Fj, Bj, J, th
S, C = sp.symbols('S C')
def canon(e):
    e = sp.expand(sp.expand_trig(e).subs({sp.tan(th): sp.sin(th)/sp.cos(th)}))
    e = sp.together(e); num, den = sp.fraction(e)
    def red(x):
        x = sp.expand(x.subs({sp.sin(th): S, sp.cos(th): C}))
        P = sp.Poly(x, C); out = 0
        for (k,), co in zip(P.monoms(), P.coeffs()): out += co*(1-S**2)**(k//2)*C**(k % 2)
        return sp.expand(out)
    return red(num), red(den)
def redSC(x):
    x = sp.expand(x); P = sp.Poly(x, C); out = 0
    for (k,), co in zip(P.monoms(), P.coeffs()): out += co*(1-S**2)**(k//2)*C**(k % 2)
    return sp.expand(out)
w = int(sys.argv[1]); t0 = time.time()
params = {6: (sp.Symbol('a'),), 2: sp.symbols('a b c0')}[w]
res = classify(w=w, params=params)
M, rhs, Bvars = res['M'], res['rhs'], res['Bvars']
Mn = sp.zeros(3,3); rn = sp.zeros(3,1)
for i in range(3):
    for j in range(3):
        n, d = canon(M[i,j]); Mn[i,j] = n/d
    n, d = canon(rhs[i]); rn[i] = n/d
detn, detd = canon(M.det()); print("det:", sp.factor(detn))
Bsol = [sp.cancel(sp.together(x)) for x in Mn.LUsolve(rn)]
subsB = dict(zip(Bvars, Bsol))
def DSC(e):
    out = sp.diff(e, S)*C - sp.diff(e, C)*S
    for k in range(J-1): out += sp.diff(e, Fj[k])*Fj[k+1]
    return out
Q = {}
for name, bs in zip(('b20','b11','b02'), Bsol):
    n, d = canon(res['subsP'][Bj[name][1]]); bp = sp.cancel((n/d).subs(subsB))
    Q[name] = sp.numer(sp.cancel(sp.together(DSC(bs) - bp)))
print("Q built", time.time()-t0, flush=True)
# structure: linear in F3?  eliminate and factor
al = {n: redSC(sp.Poly(q, Fj[3]).coeff_monomial(Fj[3])) for n, q in Q.items()}
be = {n: redSC(sp.Poly(q, Fj[3]).coeff_monomial(1)) for n, q in Q.items()}
for n in Q: print(n, "deg in F3:", sp.Poly(Q[n], Fj[3]).degree(), " alpha:", str(sp.factor(al[n]))[:150])
names = list(Q)
for n1, n2 in [(names[0], names[1]), (names[0], names[2])]:
    rel = sp.factor(redSC(al[n1]*be[n2] - al[n2]*be[n1]))
    print(f"relation {n1}-{n2}:", str(rel)[:600])
print("time", time.time()-t0)
print("\n--- Q relations directly (no f''' present):")
for n, q in Q.items():
    print(n, ":", str(sp.factor(redSC(q)))[:700])
