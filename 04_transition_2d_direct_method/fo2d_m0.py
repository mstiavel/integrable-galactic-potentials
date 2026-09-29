import sympy as sp, pickle
r = sp.Symbol('r'); alpha, beta, delta = sp.symbols('alpha beta delta')
conds, names = pickle.load(open('fo2d_m0.pkl','rb'))
B = sp.symbols('B2_0:8')
core = [sp.expand([f for f, mlt in sp.factor_list(c)[1] if f.has(B[0])][0]) for c in conds if c.has(B[0])]
print("mode-0 conditions:", len(core), "orders:", [max(k for k in range(8) if c.has(B[k])) for c in core])
# mode-2 ODE: r^2 (1+r)(1+3r) B'' + 2 r (9r^2+8r+2) B' + 2 (9r^2+4r+1) B = 0 ; substitute its solutions' jets
f = sp.Function('f')(r)
ode2 = r**2*(1+r)*(1+3*r)*sp.diff(f, r, 2) + 2*r*(9*r**2+8*r+2)*sp.diff(f, r) + 2*(9*r**2+4*r+1)*f
B2sol = sp.solve(ode2, sp.diff(f, r, 2))[0]
jets = {B[0]: f, B[1]: sp.diff(f, r)}
cur = B2sol; jets[B[2]] = cur
for k in range(3, 8):
    cur = sp.diff(cur, r).subs(sp.Derivative(f, (r, 2)), B2sol); cur = sp.simplify(cur); jets[B[k]] = cur
for c in core:
    v = sp.simplify(c.subs(jets))
    v = sp.factor(sp.numer(sp.together(v)))
    print("  mode-0 condition on mode-2 solutions:", str(v)[:300])
# test the exact 1/r^2 profile
for c in core:
    v = sp.simplify(c.subs({B[k]: sp.diff(1/r**2, r, k) for k in range(8)}))
    print("  on B=1/r^2:", sp.factor(v))
