import sympy as sp, itertools
from direct_axis import hierarchy
z, eps, s = sp.symbols('z epsilon s')
Uf = sp.Function('U')(z); Wf = sp.Function('W')(z); zd = sp.Symbol('zd')

def first_order_eqs(n, seedcoefs, A):
    rem, cs, sol = hierarchy(n, 0); free = sorted(cs.keys())
    E = sp.Rational(1,2)*zd**2 + A
    cQ = sp.expand(sum(c*E**i for i, c in enumerate(seedcoefs))*z**2)
    seed = {k: sp.Poly(cQ, zd).coeff_monomial(zd**k) for k in free}
    B1 = sp.Function('B1')(z); gam = {k: sp.Function(f'gam{k}')(z) for k in free[:-1]}
    U = A + eps*B1; W = sp.diff(U, z)/z - 3*eps*B1/z**2
    subs = {Uf: U, Wf: W}
    for k in free: subs[cs[k]] = seed[k] + (eps*gam[k] if k in gam else 0)
    eqs = []
    for m, e in rem:
        ex = sp.expand(e.subs(subs).doit())
        assert sp.simplify(ex.subs(eps, 0)) == 0
        eqs.append((m, sp.expand(sp.diff(ex, eps).subs(eps, 0))))
    return eqs, gam, B1

def leading_system(n, seedcoefs, A, shifts):
    eqs, gam, B1 = first_order_eqs(n, seedcoefs, A)
    amps = {k: sp.Symbol(f'g{k}') for k in gam}
    rep = {B1: z**s}
    for k in gam: rep[gam[k]] = amps[k]*z**(s + shifts[k])
    rows = []
    for m, e in eqs:
        ex = sp.expand(e.subs(rep).doit())
        if ex == 0: continue
        terms = sp.Add.make_args(ex)
        expo = [sp.powsimp(t).as_powers_dict().get(z, 0) for t in terms]
        # exponents are s + const; take min const
        consts = [sp.simplify(e_ - s) for e_ in expo]
        mn = min(consts)
        lead = sum(t for t, c in zip(terms, consts) if c == mn)
        lead = sp.expand(lead / z**(s + mn))
        rows.append((m, mn, lead))
    return rows, amps

if __name__ == '__main__':
    a, b, d = sp.symbols('alpha beta delta')
    for n, seed in [(4, [a, b]), (6, [a, b, d])]:
        eqs, gam, B1 = first_order_eqs(n, seed, z)
        print(f"\n=== n={n}, A=z: exponent structure (constant offsets relative to s) per equation")
        amps = {k: sp.Symbol(f'g{k}') for k in gam}
        rep = {B1: z**s}
        for k in gam: rep[gam[k]] = amps[k]*z**(s + sp.Symbol(f'p{k}'))
        for m, e in eqs:
            ex = sp.expand(e.subs(rep).doit())
            info = {}
            for t in sp.Add.make_args(ex):
                key = 'B' if not any(t.has(v) for v in amps.values()) else [str(v) for v in amps.values() if t.has(v)][0]
                c = sp.simplify(sp.powsimp(t).as_powers_dict().get(z, 0) - s)
                info.setdefault(key, set()).add(c)
            print(f"  [zdot^{m}]:", {k: sorted(v, key=lambda x: str(x)) for k, v in info.items()})
