import sympy as sp, pickle
r = sp.Symbol('r'); alpha, beta, delta = sp.symbols('alpha beta delta')
core = pickle.load(open('fo2d_m6_core.pkl','rb'))
B = sp.symbols('B6_0:8')
M, rhs = sp.linear_eq_to_matrix(core, list(B[:4]))
assert all(x == 0 for x in rhs)
from sympy.polys.matrices import DomainMatrix
K_ = sp.QQ.frac_field(r, alpha, beta, delta)
Md = DomainMatrix.from_Matrix(M).convert_to(K_)
ns = Md.nullspace().to_Matrix()
print("nullspace dimension:", ns.shape[0])
v = [sp.cancel(x) for x in ns[0, :]]
phi1 = sp.cancel(v[1]/v[0]); phi2 = sp.cancel(v[2]/v[0]); phi3 = sp.cancel(v[3]/v[0])
# consistency: (ln B)' = phi1 ; B''/B = phi1' + phi1^2 must equal phi2 ; B'''/B = phi2' + phi1*phi2 must equal phi3
c1 = sp.factor(sp.numer(sp.cancel(sp.diff(phi1, r) + phi1**2 - phi2)))
c2 = sp.factor(sp.numer(sp.cancel(sp.diff(phi2, r) + phi1*phi2 - phi3)))
print("consistency 1:", str(c1)[:600])
print("consistency 2:", str(c2)[:600])
# which (alpha,beta,delta) make them vanish identically in r?
for name, c in [('c1', c1), ('c2', c2)]:
    P = sp.Poly(sp.expand(c), r)
    sol = sp.solve(P.coeffs(), [alpha, beta, delta], dict=True)
    print(f"  {name}: seed solutions making it identically zero:", sol)
print("phi1 (= B'/B):", sp.factor(phi1))
