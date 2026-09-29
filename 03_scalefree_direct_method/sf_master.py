import sympy as sp, pickle
from sf_classify import Fj, th
S, C = sp.symbols('S C'); a, b = sp.symbols('a b')
Q = pickle.load(open('Q_w4.pkl','rb'))
def red(x):
    x = sp.expand(x); P = sp.Poly(x, C); out = 0
    for (k,), co in zip(P.monoms(), P.coeffs()): out += co*(1-S**2)**(k//2)*C**(k % 2)
    return sp.expand(out)
# master ODE: (a S^2 + b C^2)(F2 - 2 F0) - 3 (a-b) S C F1 = 0  -> F2 ; F3 by differentiation (S'=C, C'=-S)
g = a*S**2 + b*C**2
F2sol = 2*Fj[0] + 3*(a-b)*S*C*Fj[1]/g
def D(e): return sp.diff(e, S)*C - sp.diff(e, C)*S + sp.diff(e, Fj[0])*Fj[1] + sp.diff(e, Fj[1])*Fj[2]
F3sol = sp.cancel(D(F2sol).subs(Fj[2], F2sol))
print("Check: do the third-order conditions follow from the master ODE?")
for name, q in Q.items():
    v = sp.cancel(sp.together(q.subs(Fj[3], F3sol).subs(Fj[2], F2sol)))
    v = red(sp.numer(v))
    print("  ", name, ":", "0" if v == 0 else "NONZERO")
# solve the master ODE exactly
f = sp.Function('f')(th); ab = sp.symbols('kappa')   # kappa = b/a
ode = sp.Eq((sp.sin(th)**2 + ab*sp.cos(th)**2)*(f.diff(th,2) - 2*f) - 3*(1-ab)*sp.sin(th)*sp.cos(th)*f.diff(th), 0)
try:
    sol = sp.dsolve(ode, f)
    print("general solution:", sol)
except Exception as e:
    print("dsolve failed:", e)
