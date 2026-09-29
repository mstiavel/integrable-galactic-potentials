import sympy as sp, pickle
from sf_classify import Fj, th
S, C = sp.symbols('S C'); a, b = sp.symbols('a b')
Q = pickle.load(open('Q_w4.pkl','rb'))
def red(x):
    x = sp.expand(x); P = sp.Poly(x, C); out = 0
    for (k,), co in zip(P.monoms(), P.coeffs()): out += co*(1-S**2)**(k//2)*C**(k % 2)
    return sp.expand(out)
g = a*S**2 + b*C**2
F2sol = 2*Fj[0] + 3*(a-b)*S*C*Fj[1]/g
def D(e): return sp.diff(e, S)*C - sp.diff(e, C)*S + sp.diff(e, Fj[0])*Fj[1] + sp.diff(e, Fj[1])*Fj[2]
F3sol = sp.cancel(D(F2sol).subs(Fj[2], F2sol))
rels = {}
for name, q in Q.items():
    v = sp.cancel(sp.together(q.subs(Fj[3], F3sol).subs(Fj[2], F2sol)))
    v = sp.factor(red(sp.numer(v)))
    rels[name] = v
    print(f"{name}: {v}\n")
pickle.dump(rels, open('rels_master.pkl','wb'))
