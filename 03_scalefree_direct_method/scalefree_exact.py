import sympy as sp, itertools, sys
from direct_axis import hierarchy
z, w, t = sp.symbols('z w t')
Uf = sp.Function('U')(z); Wf = sp.Function('W')(z)

def exact_conditions(n, orbit):
    rem, cs, sol = hierarchy(n, 0)
    free = sorted(cs.keys())
    Cs = {k: sp.Symbol(f'C{k}') for k in free}
    # scale-free mu=1 potential r(1 + t P2): axis U=(1+t)z, W=(1-2t)/z ; equator V=(1-t/2)x, W=(1+5t/2)/x
    if orbit == 'axis': U = (1+t)*z; W = (1-2*t)/z
    else:               U = (1-t/2)*z; W = (1+sp.Rational(5,2)*t)/z
    subs = {Uf: U, Wf: W}
    for k in free: subs[cs[k]] = Cs[k]*z**(w - 1 - sp.Rational(k, 2))
    conds = []
    for m, e in rem:
        ex = sp.simplify(e.subs(subs).doit())
        ex = sp.expand(sp.powsimp(sp.expand(ex), force=True))
        # factor out the common power of z
        ex = sp.simplify(ex / z**(w - 4))  # placeholder normalisation
        ex = sp.expand(sp.powsimp(ex, force=True))
        # collect coefficient by removing any remaining z power
        pows = {sp.powsimp(term).as_powers_dict().get(z, 0) for term in sp.Add.make_args(ex)}
        assert len(pows) <= 1, (m, pows)
        if ex != 0:
            conds.append(sp.factor(sp.expand(ex / z**list(pows)[0])))
    return conds, Cs

for n in [4, 6, 8]:
    print(f"\n======== n = {n}")
    res = {}
    for orbit in ['axis', 'eq']:
        conds, Cs = exact_conditions(n, orbit)
        # solve: nontrivial (C_k not all zero) solutions; treat as polynomial system in C's, w, t
        sols = sp.solve(conds, list(Cs.values()) + [w], dict=True)
        tvals = set()
        print(f"  {orbit}: {len(conds)} conditions; solution branches:")
        for s_ in sols:
            Cvals = {k: s_.get(Cs[k], Cs[k]) for k in Cs}
            if all(v == 0 for v in Cvals.values()): continue
            print("     ", {str(k): v for k, v in s_.items()})
        # additionally solve for t explicitly with the top C nonzero
        top = Cs[max(Cs)]
        sols_t = sp.solve(conds + [top - 1], list(Cs.values()) + [w, t], dict=True)
        ts = sorted({sp.nsimplify(s_[t]) for s_ in sols_t if t in s_})
        print(f"     t values with top coefficient nonzero: {ts}")
        res[orbit] = ts
    sys.stdout.flush()
