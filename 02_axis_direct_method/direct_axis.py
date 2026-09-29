import sympy as sp

z, zd, E, kappa, C = sp.symbols('z zdot E kappa C')
U = sp.Function('U')(z)
W = sp.Function('W')(z)

def D(f):
    """time derivative along the axial orbit: D = zdot d/dz - U' d/dzdot"""
    return zd*sp.diff(f, z) - sp.diff(U, z)*sp.diff(f, zd)

def hierarchy(n, parity=None):
    """Return (constraints, coefficient functions) for a conserved quadratic form
       Q = a xi^2 + 2 b xi xidot + c xidot^2 with deg_zdot (a,b,c) = (n, n-1, n-2)."""
    cs = {k: sp.Function(f'c{k}')(z) for k in range(0, n-1)}
    bs = {k: sp.Function(f'b{k}')(z) for k in range(0, n)}
    as_ = {k: sp.Function(f'a{k}')(z) for k in range(0, n+1)}
    if parity is not None:   # keep only zdot^k with k = parity (mod 2) for c, opposite for b, parity for a
        cs = {k: v for k, v in cs.items() if k % 2 == parity}
        bs = {k: v for k, v in bs.items() if k % 2 != parity}
        as_ = {k: v for k, v in as_.items() if k % 2 == parity}
    c = sum(v*zd**k for k, v in cs.items())
    b = sum(v*zd**k for k, v in bs.items())
    a = sum(v*zd**k for k, v in as_.items())
    eqs = [sp.expand(D(c) + 2*b), sp.expand(D(b) + a - c*W), sp.expand(D(a) - 2*b*W)]
    # solve the first two for b_k, a_k in terms of c_k
    sol = {}
    e1 = sp.Poly(eqs[0], zd)
    for k, bk in bs.items():
        coeff = e1.coeff_monomial(zd**k)
        sol[bk] = sp.solve(coeff, bk)[0]
    e2 = sp.Poly(sp.expand(eqs[1].subs(sol)), zd)
    for k, ak in as_.items():
        coeff = e2.coeff_monomial(zd**k)
        sol[ak] = sp.solve(coeff, ak)[0]
    # remaining constraints from third equation and leftover monomials of the first two
    rem = []
    for e in eqs:
        p = sp.Poly(sp.expand(e.subs(sol).doit()), zd)
        for m in p.monoms():
            cf = sp.simplify(p.coeff_monomial(zd**m[0]))
            if cf != 0: rem.append((m[0], cf))
    return rem, cs, sol

if __name__ == '__main__':
    for n, par in [(2, 0), (3, 1), (4, 0)]:
        rem, cs, sol = hierarchy(n, par)
        print(f"\n===== n = {n}: free functions {list(cs.values())}")
        for m, cf in rem:
            print(f"  [zdot^{m}]  ", sp.factor(cf), " = 0")
