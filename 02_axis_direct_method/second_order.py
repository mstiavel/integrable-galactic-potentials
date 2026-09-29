import sympy as sp, itertools, sys, time
from direct_axis import hierarchy

z, eps = sp.symbols('z epsilon')
Uf = sp.Function('U')(z); Wf = sp.Function('W')(z)
rem, cs, sol = hierarchy(4, 0); c0f, c2f = cs[0], cs[2]
E3 = rem[1][1]; E1 = rem[2][1]

def run(N, M, alpha, beta, verbose=True, n2=False):
    t0 = time.time()
    # cored isothermal on axis: A = 1/2 ln(1+z^2), Taylor to high order
    A = sp.series(sp.log(1+z**2)/2, z, 0, N+8).removeO()
    # unknown series
    b1 = sp.symbols(f'b1_1:{N+1}'); b2 = sp.symbols(f'b2_1:{N+1}')
    B1 = sum(c*z**(i+1) for i, c in enumerate(b1)); B2 = sum(c*z**(i+1) for i, c in enumerate(b2))
    gam = {}
    for tag in ['1ax', '1eq', '2ax', '2eq']:
        gam[tag] = sp.symbols(f'g{tag}_m2 g{tag}_m1') + sp.symbols(f'g{tag}_0:{N+1}')
    def gser(tag): return sum(c*z**(i-2) for i, c in enumerate(gam[tag]))
    ks = {tag: sp.symbols(f'k{tag}_0:3') for tag in ['1ax', '1eq', '2ax', '2eq']}
    def kpoly(tag): return 0 if n2 else ks[tag][0] + ks[tag][1]*z + ks[tag][2]*z**2

    eqs = {1: [], 2: []}
    for orbit in ['ax', 'eq']:
        Bfull = eps*B1 + eps**2*B2
        if orbit == 'ax':
            U = A + Bfull; W = sp.diff(U, z)/z - 3*Bfull/z**2
        else:
            U = A - Bfull/2; W = sp.diff(U, z)/z + 3*Bfull/z**2
        c0 = (alpha + beta*A)*z**2 + eps*gser('1'+orbit) + eps**2*gser('2'+orbit)
        c2 = beta*z**2/2 + eps*kpoly('1'+orbit) + eps**2*kpoly('2'+orbit)
        subs = {c0f: c0, c2f: c2, Uf: U, Wf: W}
        for E in (E3, E1):
            ex = E.subs(subs).doit()
            ex = sp.expand(ex * z**6)          # clear the explicit 1/z powers
            ex = sp.expand(ex)
            # expand in eps to order 2
            pe = sp.Poly(ex, eps)
            for order in (1, 2):
                co = pe.coeff_monomial(eps**order)
                co = sp.expand(co)
                pz = sp.Poly(co, z)
                for mon in pz.monoms():
                    if mon[0] <= M + 6:
                        eqs[order].append(pz.coeff_monomial(z**mon[0]))
    if verbose: print(f"  built equations in {time.time()-t0:.0f}s: O(eps) {len(eqs[1])}, O(eps^2) {len(eqs[2])}")
    # ---- first order: linear homogeneous in x1
    x1 = list(b1) + list(gam['1ax']) + list(gam['1eq']) + ([] if n2 else list(ks['1ax']) + list(ks['1eq']))
    M1, r1 = sp.linear_eq_to_matrix(eqs[1], x1)
    assert all(r == 0 for r in r1)
    ns1 = M1.nullspace()
    if verbose: print(f"  first order: {len(ns1)} free directions")
    tpar = sp.symbols(f't0:{len(ns1)}')
    gen1 = sum((tp*v for tp, v in zip(tpar, ns1)), sp.zeros(len(x1), 1))
    sub1 = {xi: gen1[i] for i, xi in enumerate(x1)}
    # which directions involve low-order b1 coefficients (genuine, not truncation artifacts)?
    lowb = [sub1[b] for b in b1[:M]]
    if verbose: print("  low-order B1 coefficients in terms of free params:", lowb)
    # ---- second order: linear in x2 with sources quadratic in t
    x2 = list(b2) + list(gam['2ax']) + list(gam['2eq']) + ([] if n2 else list(ks['2ax']) + list(ks['2eq']))
    eqs2 = [sp.expand(e.subs(sub1)) for e in eqs[2]]
    M2, r2 = sp.linear_eq_to_matrix(eqs2, x2)
    # solvability: left nullspace of M2 dotted with r2 must vanish
    LN = M2.T.nullspace()
    conds = [sp.expand((n.T*r2)[0]) for n in LN]
    conds = [c for c in conds if c != 0]
    if verbose: print(f"  second order: {len(LN)} left-null directions, {len(conds)} nontrivial compatibility conditions")
    return conds, tpar, lowb, sub1, b1


def b1_family_dim(conds, tpar, lowb, b1, M):
    """dimension of the space of low-order B1 coefficient vectors compatible with the conditions (quadratic system)."""
    sols = sp.solve(conds, tpar, dict=True) if conds else [{}]
    dims = []
    for s_ in sols:
        vec = sp.Matrix([sp.simplify(x.subs(s_)) for x in lowb])
        free = sorted(vec.free_symbols, key=str)
        J = vec.jacobian(free) if free else sp.zeros(len(vec), 0)
        # generic rank (substitute random values)
        import random
        Jn = J.subs({f: sp.Rational(random.randint(2, 50), random.randint(2, 50)) for f in free})
        dims.append((Jn.rank(), [sp.factor(v) for v in vec]))
    return dims

if __name__ == '__main__':
    N, M = int(sys.argv[1]), int(sys.argv[2])
    for alpha, beta, n2 in [(1, 0, True), (1, 0, False), (1, sp.Rational(1,3), False)]:
        print(f"\n=== alpha={alpha}, beta={beta}, n2-restricted={n2}, N={N}, M={M}")
        conds, tpar, lowb, sub1, b1 = run(N, M, alpha, beta, n2=n2)
        for d, vec in b1_family_dim(conds, tpar, lowb, b1, M):
            print(f"  surviving B1 family: dimension {d}; coefficients b1..b{M}:", vec)
        sys.stdout.flush()
