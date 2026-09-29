import sympy as sp, sys, itertools
from sf_direct import system, th, f
S, C = sp.symbols('S C')   # sin theta, cos theta

def canon(expr):
    """expression in sin(th), cos(th) -> polynomial in S, C with C^2 reduced by 1 - S^2"""
    e = sp.expand(expr.subs({sp.sin(th): S, sp.cos(th): C}))
    P = sp.Poly(e, C)
    out = 0
    for (k,), co in zip(P.monoms(), P.coeffs()):
        out += co * (1 - S**2)**(k//2) * C**(k % 2)
    return sp.expand(out)

def solve_fixed(w, tval, D=8, verbose=True):
    nparams = {2: 3, 4: 2, 6: 1}[w]
    params = sp.symbols(f'a0:{nparams}')
    eqs, (b20, b11, b02, c) = system(w, params)
    fexpr = 1 + tval*(3*sp.cos(th)**2 - 1)/2
    # trig-polynomial ansatz for the unknown functions
    coeffs = []
    def ansatz(name):
        terms = 0
        for i in range(D+1):
            for j in range(0, 2):              # C^0 or C^1 (C^2 reduced)
                if i + j > D: continue
                sym = sp.Symbol(f'{name}_{i}_{j}'); coeffs.append(sym)
                terms += sym*sp.sin(th)**i*sp.cos(th)**j
        return terms
    A20, A11, A02, Ac = ansatz('b20'), ansatz('b11'), ansatz('b02'), ansatz('c')
    lin_eqs = []
    for key, e in eqs:
        ex = e.subs(f, fexpr).subs({b20: A20, b11: A11, b02: A02, c: Ac}).doit()
        ex = canon(ex)
        P = sp.Poly(ex, S, C)
        lin_eqs += [co for co in P.coeffs()]
    unknowns = coeffs + list(params)
    M, rhs = sp.linear_eq_to_matrix(lin_eqs, unknowns)
    assert all(r == 0 for r in rhs)
    ns = M.nullspace()
    # solutions with nonzero Killing part
    good = []
    for v in ns:
        kpart = v[len(coeffs):, 0]
        if any(x != 0 for x in kpart): good.append(v)
    # dimension of Killing-parameter projection
    if good:
        Kmat = sp.Matrix.hstack(*[v[len(coeffs):, 0] for v in good])
        kdim = Kmat.rank()
    else:
        kdim = 0
    if verbose:
        print(f"  w={w}, t={tval}: nullspace dim {len(ns)}, solutions with nonzero Killing tensor: {len(good)}, rank of Killing projection: {kdim}")
        for v in good[:4]:
            kp = {str(p): v[len(coeffs)+i, 0] for i, p in enumerate(params)}
            print("     Killing params:", kp)
    return kdim, good, params, coeffs, (A20, A11, A02, Ac)

if __name__ == '__main__':
    tvals = [sp.Integer(0), sp.Rational(-2, 5), sp.Rational(1, 2), sp.Rational(3, 10)]
    for w in [4, 6, 2]:
        print(f"\n===== weight {w}")
        for tv in tvals:
            solve_fixed(w, tv)
            sys.stdout.flush()
