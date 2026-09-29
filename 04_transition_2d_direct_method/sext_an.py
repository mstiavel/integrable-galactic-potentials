import sympy as sp, pickle, sys
r = sp.Symbol('r')
m = int(sys.argv[1]); l = m
conds, params = pickle.load(open(f'fo2d_hernquist_n6_m{m}.pkl','rb'))
vals = {p: v for p, v in zip(params, [sp.Rational(3,7), sp.Rational(2,5), sp.Rational(5,11), sp.Rational(7,13), sp.Rational(4,9), sp.Rational(6,17)])}
B = sp.symbols(f'B{l}_0:9')
core = []
for c in conds:
    c = sp.expand(c.subs(vals))
    fs = [f for f, mlt in sp.factor_list(c)[1] if f.has(B[0])]
    if fs: core.append(sp.expand(fs[0]))
kmax = max(max(k for k in range(9) if c.has(B[k])) for c in core)
print(f"mode {m}: {len(core)} conditions, max order {kmax}")
M, rhs = sp.linear_eq_to_matrix(core, list(B[:kmax+1]))
from sympy.polys.matrices import DomainMatrix
K_ = sp.QQ.frac_field(r)
ns = DomainMatrix.from_Matrix(M).convert_to(K_).nullspace().to_Matrix()
dim = ns.shape[0]; print("  jet-nullspace dimension:", dim)
if dim == 1:
    v = [sp.cancel(x) for x in ns[0, :]]
    phi = [sp.cancel(v[k]/v[0]) for k in range(kmax+1)]
    print("  B'/B =", sp.factor(phi[1]))
    cur = phi[1]
    for k in range(2, kmax+1):
        nxt = sp.cancel(sp.diff(cur, r) + phi[1]*cur)
        print(f"  consistency at order {k}:", sp.factor(sp.numer(sp.cancel(nxt - phi[k]))))
        cur = phi[k]
elif dim == 2:
    Kb = sp.QQ.frac_field(r, *B[:2])
    sub = M[:, 2:kmax+1]; rr = -(M[:, 0:2]*sp.Matrix(B[:2]))
    X = DomainMatrix.from_Matrix(sub).convert_to(Kb).lu_solve(DomainMatrix.from_Matrix(rr).convert_to(Kb)).to_Matrix()
    B2e = sp.cancel(X[0]); a = sp.cancel(B2e.subs(B[1], 0)/B[0]); b = sp.cancel(B2e.subs(B[0], 0)/B[1])
    print("  B'' = a B + b B' with a =", sp.factor(a), "; b =", sp.factor(b))
    # consistency of the third jet
    B3e = sp.cancel(X[1]); B3c = sp.expand(sp.diff(a, r)*B[0] + a*B[1] + sp.diff(b, r)*B[1] + b*B2e)
    print("  consistency:", sp.factor(sp.numer(sp.cancel(B3c - B3e))))
    s = sp.Symbol('s'); print("  indicial at r=0:", sp.factor(sp.limit(sp.expand(s*(s-1) - b*r*s - a*r**2), r, 0)))
else:
    print("  larger solution space")
