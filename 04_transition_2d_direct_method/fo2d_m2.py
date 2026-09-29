import sympy as sp, pickle
r = sp.Symbol('r'); alpha, beta, delta = sp.symbols('alpha beta delta')
conds, names = pickle.load(open('fo2d_m2.pkl','rb'))
B = sp.symbols('B2_0:8')
core = [sp.expand([f for f, mlt in sp.factor_list(c)[1] if f.has(B[0])][0]) for c in conds]
M, rhs = sp.linear_eq_to_matrix(core, list(B[:4]))
from sympy.polys.matrices import DomainMatrix
K_ = sp.QQ.frac_field(r, alpha, beta, delta, *B[:2])
Md = DomainMatrix.from_Matrix(M).convert_to(K_)
# solve for B2, B3 in terms of B0, B1
sub = M[:, 2:4]; rr = -(M[:, 0:2]*sp.Matrix(B[:2]))
X = DomainMatrix.from_Matrix(sub).convert_to(K_).lu_solve(DomainMatrix.from_Matrix(rr).convert_to(K_)).to_Matrix()
B2e = sp.cancel(X[0]); B3e = sp.cancel(X[1])
a = sp.cancel(B2e.subs(B[1], 0)/B[0]); b = sp.cancel(B2e.subs(B[0], 0)/B[1])
print("B'' = a B + b B' with a =", sp.factor(a), " b =", sp.factor(b))
B3c = sp.expand(sp.diff(a, r)*B[0] + a*B[1] + sp.diff(b, r)*B[1] + b*B2e)
cons = sp.factor(sp.numer(sp.cancel(B3c - B3e)))
print("consistency (should vanish for a 2-dim solution space):", str(cons)[:400])
f = sp.Function('f')(r)
ode = sp.Eq(sp.diff(f, r, 2), a*f + b*sp.diff(f, r))
print("second-order ODE:", sp.factor(sp.diff(f,r,2) - a*f - b*sp.diff(f,r)))
try:
    print("solution:", sp.dsolve(ode, f))
except Exception as e:
    print("dsolve failed:", e)
# exponents at r=0
s = sp.Symbol('s')
ind = sp.factor(sp.limit(sp.expand((s*(s-1) - b*r*s - a*r**2)), r, 0))
print("indicial at r=0:", ind)
