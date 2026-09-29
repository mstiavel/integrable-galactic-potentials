import sympy as sp, sys, time
from so2d_step2 import load_I1, S2_modes, r, pr, pt, alpha, beta, delta, b, A, Ap, I0
from fo2d import legendre_modes
m = int(sys.argv[1]); N = int(sys.argv[2]); t0 = time.time()
quadratic_seed = len(sys.argv) > 3 and sys.argv[3] == 'quad'
Q2, Q0, Qm2 = load_I1(); S2 = S2_modes(Q2, Q0, Qm2)[m]
# unknown quartic Q_m: c_{jk} = r^{-k} * sum_{n=0}^{N} a_n r^n ; second-order potential C_l regular: sum_{n=2}^{N} g_n r^n
unk = []; Q = 0
for j in range(5):
    for k in range(5-j):
        if (j+k) % 2 == 0:
            ser = 0
            for n in range(N+1):
                s = sp.Symbol(f'a{j}{k}_{n}'); unk.append(s); ser += s*r**n
            Q += ser*r**(-k)*pr**j*pt**k
Cl = {}
for l in (0, 2, 4):
    ser = 0
    for n in range(2, N+1):
        s = sp.Symbol(f'g{l}_{n}'); unk.append(s); ser += s*r**n
    Cl[l] = ser
H2m = sum(Cl[l]*legendre_modes(l).get(m, 0) for l in Cl)
lhs = pr*sp.diff(Q, r) + (pt/r**2)*(sp.I*m)*Q - (Ap - pt**2/r**3)*sp.diff(Q, pr)
srcH2 = -sp.diff(I0, pr)*sp.diff(H2m, r) - sp.diff(I0, pt)*(sp.I*m*H2m)
eq = lhs + S2 + srcH2
eq = sp.series(sp.expand(eq), r, 0, N-2).removeO()   # truncate consistently
eq = sp.expand(eq)
P = sp.Poly(eq, pr, pt)
lin = []
for co in P.coeffs():
    co = sp.expand(co); Pr = sp.Poly(co, r, 1/r) if False else None
    # collect powers of r (Laurent): multiply by large power
    cs = sp.expand(co*r**12); Pp = sp.Poly(cs, r)
    lin += list(Pp.coeffs())
lin = [e for e in lin if e != 0]
M, rhs = sp.linear_eq_to_matrix(lin, unk)
print(f"mode {m}, N={N}: {len(lin)} equations, {len(unk)} unknowns, built in {time.time()-t0:.0f}s", flush=True)
# solvability: rank[M|rhs] == rank[M] ?  (rhs is in b^2 alpha,beta,delta: check generically by random numeric values)
import random
vals = {alpha: sp.Rational(3,7), beta: 0 if quadratic_seed else sp.Rational(2,5), delta: 0 if quadratic_seed else sp.Rational(5,11), b: 1}
Mn = M.subs(vals); rn = rhs.subs(vals)
rk = Mn.rank(); rka = Mn.row_join(rn).rank()
print(f"  rank M = {rk}, rank [M|rhs] = {rka}  ->", "SOLVABLE (regular second-order solution exists to this order)" if rk == rka else "NOT SOLVABLE (second-order obstruction)", flush=True)
# also with source switched off (control: must be solvable)
print("  control without source:", Mn.rank() == Mn.row_join(sp.zeros(Mn.shape[0],1)).rank())
# freedom in the second-order potential: dimension of the projection of the solution space onto the g-coefficients (low orders only, n <= N-4)
gcols = [i for i, u in enumerate(unk) if str(u).startswith('g') and int(str(u).split('_')[1]) <= N-4]
ns = Mn.nullspace()
G = sp.Matrix.hstack(*[v.extract(gcols, [0]) for v in ns]) if ns else sp.zeros(len(gcols), 0)
print(f"  free directions in low-order g's (potential freedom, homogeneous): {G.rank()} of {len(gcols)}")
if rk == rka:
    sol = Mn.gauss_jordan_solve(rn)[0]
    gpart = {str(unk[i]): sp.nsimplify(sol[i]) for i in gcols}
    print("  particular solution g's:", {k: v for k, v in gpart.items() if v != 0 and not v.free_symbols})
