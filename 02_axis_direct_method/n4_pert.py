import sympy as sp, pickle
from direct_axis import hierarchy
z, eps = sp.symbols('z epsilon')
alpha, beta = sp.symbols('alpha beta')
Uf = sp.Function('U')(z); Wf = sp.Function('W')(z)

rem, cs, sol = hierarchy(4, 0)
c0f, c2f = cs[0], cs[2]
E3 = rem[1][1]   # [zdot^3]
E1 = rem[2][1]   # [zdot^1]
assert rem[0][0] == 5 and rem[1][0] == 3 and rem[2][0] == 1

A = -1/(1+z)                      # Hernquist on axis (cusp at 0, spherical)
B1 = sp.Function('B1')(z)

def orbit_setup(orbit):
    g = sp.Function('gamma_'+orbit)(z)
    k0, k1, k2 = sp.symbols(f'k0_{orbit} k1_{orbit} k2_{orbit}')
    if orbit == 'axis':
        U = A + eps*B1; W = sp.diff(U, z)/z - 3*eps*B1/z**2
    else:
        U = A - eps*B1/2; W = sp.diff(U, z)/z + 3*eps*B1/z**2
    U0 = A
    c0 = (alpha + beta*U0)*z**2 + eps*g
    c2 = beta*z**2/2 + eps*(k0 + k1*z + k2*z**2)
    subs = {c0f: c0, c2f: c2, Uf: U, Wf: W}
    e3 = sp.expand(E3.subs(subs).doit()); e1 = sp.expand(E1.subs(subs).doit())
    e3_0 = sp.simplify(e3.subs(eps, 0)); e1_0 = sp.simplify(e1.subs(eps, 0))
    e3_1 = sp.simplify(sp.diff(e3, eps).subs(eps, 0)); e1_1 = sp.simplify(sp.diff(e1, eps).subs(eps, 0))
    return g, (k0, k1, k2), (e3_0, e1_0), (e3_1, e1_1)

out = {}
for orbit in ['axis', 'eq']:
    g, ks, zero, first = orbit_setup(orbit)
    print(f"\n=== {orbit}: zeroth-order residuals (must vanish):", zero)
    e3_1, e1_1 = first
    # e3_1 is linear in gamma''' and B1 (up to B1'''); solve e3_1 for gamma''' ... instead eliminate gamma:
    # both are linear in gamma and derivatives; use e3_1 to express gamma''' and substitute into d^2/dz^2 ... simpler:
    # treat as linear system: solve e1_1 for gamma'' ? Let's print orders.
    print("  order in gamma:", sp.ode_order(e3_1, g), sp.ode_order(e1_1, g), " order in B1:", sp.ode_order(e3_1, B1), sp.ode_order(e1_1, B1))
    out[orbit] = (g, ks, e3_1, e1_1)
pickle.dump(out, open('n4_first_order.pkl', 'wb'))
