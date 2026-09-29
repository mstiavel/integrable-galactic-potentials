import sympy as sp, pickle, time
from sf_classify import Fj
S, C = sp.symbols('S C'); a, b = sp.symbols('a b')
Q = pickle.load(open('Q_w4.pkl','rb'))
t0 = time.time()
def red(x):
    x = sp.expand(x); P = sp.Poly(x, C); out = 0
    for (k,), co in zip(P.monoms(), P.coeffs()): out += co*(1-S**2)**(k//2)*C**(k % 2)
    return sp.expand(out)
al, be = {}, {}
for name, q in Q.items():
    P = sp.Poly(q, Fj[3])
    print(name, "degree in f''':", P.degree(), " monomials:", P.monoms())
    al[name] = red(P.coeff_monomial(Fj[3])); be[name] = red(P.coeff_monomial(1))
    print("   alpha factor:", sp.factor(al[name])[:0] if False else str(sp.factor(al[name]))[:160])
# eliminate f''': two algebraic relations
names = list(Q)
Rel = {}
for n1, n2 in [(names[0], names[1]), (names[0], names[2])]:
    rel = red(al[n1]*be[n2] - al[n2]*be[n1])
    rel = sp.factor(rel)
    Rel[(n1, n2)] = rel
    P = sp.Poly(sp.expand(rel), Fj[0], Fj[1], Fj[2])
    print(f"relation {n1}-{n2}: total degree in (f,f',f'') {P.total_degree()}, #terms {len(P.terms())}, factors (short):", str(rel)[:300])
pickle.dump((al, be, Rel), open('rel_w4.pkl','wb'))
print("time", time.time()-t0)
