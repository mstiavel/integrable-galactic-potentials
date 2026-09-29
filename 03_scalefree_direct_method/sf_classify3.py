import sympy as sp, pickle, time
from sf_classify import *
from sf_classify2 import canon, S, C, DSC, Bsol, res, subsB
t0 = time.time()
Bvars = res['Bvars']
# b' expressions in canonical (S,C) form with b's substituted by their solutions
bprime = {}
for name in ('b20','b11','b02'):
    n, d = canon(res['subsP'][Bj[name][1]]); bprime[name] = sp.cancel((n/d).subs(subsB))
Q = {}
for name, bs in zip(('b20','b11','b02'), Bsol):
    Q[name] = sp.numer(sp.cancel(sp.together(DSC(bs) - bprime[name])))
    print(name, "Q order in f:", max(k for k in range(J) if Q[name].has(Fj[k])), "size", len(str(Q[name])), flush=True)
pickle.dump(Q, open('Q_w4.pkl','wb'))
a, b = res['params']
def red(x):
    x = sp.expand(x); P = sp.Poly(x, C); out = 0
    for (k,), co in zip(P.monoms(), P.coeffs()): out += co*(1-S**2)**(k//2)*C**(k % 2)
    return sp.expand(out)
def check(fexpr, av, bv, label):
    fj = {Fj[k]: sp.diff(fexpr, th, k).subs({sp.sin(th): S, sp.cos(th): C}) for k in range(J)}
    out = []
    for name in Q:
        v = red(Q[name].subs({a: av, b: bv}).subs(fj))
        out.append('0' if v == 0 else 'X')
    print(f"  {label:34s} Q(b20,b11,b02): {out}")
check(sp.Integer(1), 1, 1, 'spherical a=b')
check(sp.Integer(1), 1, 0, 'spherical a!=b (must fail)')
check(1 + sp.Rational(1,2)*(3*sp.cos(th)**2-1)/2, 1, 0, 't=1/2, (1,0)')
check(1 + sp.Rational(1,2)*(3*sp.cos(th)**2-1)/2, 0, 1, 't=1/2, (0,1) (must fail)')
check(1 - sp.Rational(2,5)*(3*sp.cos(th)**2-1)/2, 0, 1, 't=-2/5, (0,1)')
check(1 - sp.Rational(2,5)*(3*sp.cos(th)**2-1)/2, 1, 1, 't=-2/5, (1,1) (must fail)')
check(1 + sp.Rational(3,10)*(3*sp.cos(th)**2-1)/2, 1, 0, 't=0.3 (must fail)')
print("time", time.time()-t0)
