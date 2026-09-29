import sympy as sp, pickle, sys
r = sp.Symbol('r'); alpha, beta, delta = sp.symbols('alpha beta delta')
m = int(sys.argv[1]); l = int(sys.argv[2])
conds, names = pickle.load(open(f'fo2d_m{m}.pkl','rb'))
B = sp.symbols(f'B{l}_0:8')
core = [sp.expand([f for f, mlt in sp.factor_list(c)[1] if f.has(B[0])][0]) for c in conds]
orders = [max(k for k in range(8) if c.has(B[k])) for c in core]
print("orders:", orders)
kmax = max(orders)
M, rhs = sp.linear_eq_to_matrix(core, list(B[:kmax+1]))
from sympy.polys.matrices import DomainMatrix
K_ = sp.QQ.frac_field(r, alpha, beta, delta)
ns = DomainMatrix.from_Matrix(M).convert_to(K_).nullspace().to_Matrix()
print("nullspace dimension (jet space):", ns.shape[0])
if ns.shape[0] == 1:
    v = [sp.cancel(x) for x in ns[0, :]]
    phi = [sp.cancel(v[k]/v[0]) for k in range(kmax+1)]
    print("B'/B =", sp.factor(phi[1]))
    # consistency
    ok = True
    cur = phi[1]
    for k in range(2, kmax+1):
        nxt = sp.cancel(sp.diff(cur, r) + phi[1]*cur)
        c = sp.factor(sp.numer(sp.cancel(nxt - phi[k])))
        print(f"  consistency at order {k}:", 0 if c == 0 else str(c)[:300])
        cur = phi[k]
else:
    print("larger solution space: analyse further")
