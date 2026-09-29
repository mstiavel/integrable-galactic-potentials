import sympy as sp, pickle, time, sys
from sf_classify import *
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

t0 = time.time()
res = classify()
M, rhs, Bvars = res['M'], res['rhs'], res['Bvars']
# canonical forms
Mc = M.applyfunc(lambda e: sp.cancel(sp.Mul(*[sp.Rational(1)])*e))
Mn = sp.zeros(3,3); Dn = sp.zeros(3,3)
for i in range(3):
    for j in range(3):
        n, d = canon(M[i,j]); Mn[i,j] = n/d
rn = sp.zeros(3,1)
for i in range(3):
    n, d = canon(rhs[i]); rn[i] = n/d
Bsol = Mn.LUsolve(rn)
Bsol = [sp.cancel(sp.together(x)) for x in Bsol]
print(f"b's solved in {time.time()-t0:.0f}s; sizes", [len(str(x)) for x in Bsol], flush=True)
pickle.dump(Bsol, open('Bsol_w4.pkl','wb'))
# derivative rules in (S,C) variables: d/dθ S = C, d/dθ C = -S
def DSC(e):
    out = sp.diff(e, S)*C - sp.diff(e, C)*S
    for k in range(J-1): out += sp.diff(e, Fj[k])*Fj[k+1]
    return out
# remaining relations: R3' and R1'' must vanish after substituting b's (their derivatives via subsP, then b's)
subsB = dict(zip(Bvars, Bsol))
def full(e):
    # e in jets: replace b' by subsP (which are in theta form) -> canon -> substitute b's
    e2 = sp.expand(e.subs(res['subsP']))
    n, d = canon(e2)
    return sp.cancel(sp.together((n/d).subs(subsB)))
R3, R1p = res['R3'], res['R1p']
# D(R3): first canonicalize R3 with b's still symbolic, then differentiate in (S,C, jets incl. b jets)
def Dfull(e):
    n, d = canon(sp.expand(e.subs(res['subsP'])))   # now in S,C, F, b0's
    ex = n/d
    out = sp.diff(ex, S)*C - sp.diff(ex, C)*S
    for k in range(J-1): out += sp.diff(ex, Fj[k])*Fj[k+1]
    for name in ('b20','b11','b02'):
        bp = res['subsP'][Bj[name][1]]; nb, db = canon(bp)
        out += sp.diff(ex, Bj[name][0])*nb/db
    return sp.cancel(sp.together(out.subs(subsB)))
ODE_A = sp.numer(Dfull(R3)); ODE_B = sp.numer(Dfull(R1p))
print("ODE_A order in f:", max(k for k in range(J) if ODE_A.has(Fj[k])), " size", len(str(ODE_A)))
print("ODE_B order in f:", max(k for k in range(J) if ODE_B.has(Fj[k])), " size", len(str(ODE_B)))
pickle.dump((ODE_A, ODE_B), open('fODEs_w4.pkl','wb'))
# checks: spherical f=1 with a=b, parabolic f=1+C^2 with (a,b)=(1,0), f = 1 - (3/2) ... t=-2/5 profile with (0,1)
a, b = res['params']
def check(fexpr, av, bv, label):
    fj = {Fj[k]: sp.diff(fexpr, th, k).subs({sp.sin(th): S, sp.cos(th): C}) for k in range(J)}
    for name, ode in (('A', ODE_A), ('B', ODE_B)):
        v = sp.expand(ode.subs({a: av, b: bv}).subs(fj))
        v = sp.simplify(v.subs(C**2, 1 - S**2))
        print(f"  check {label} ODE_{name}:", sp.simplify(v) == 0 or sp.simplify(v))
check(sp.Integer(1), 1, 1, 'spherical a=b')
check(sp.Integer(1), 1, 0, 'spherical a!=b (must fail)')
check(1 + sp.Rational(1,2)*(3*sp.cos(th)**2-1)/2, 1, 0, 't=1/2, (1,0)')
check(1 - sp.Rational(2,5)*(3*sp.cos(th)**2-1)/2, 0, 1, 't=-2/5, (0,1)')
print("total", time.time()-t0)
