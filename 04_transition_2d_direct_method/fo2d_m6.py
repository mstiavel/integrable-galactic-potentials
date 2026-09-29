import sympy as sp, pickle
r = sp.Symbol('r'); alpha, beta, delta = sp.symbols('alpha beta delta')
conds, names = pickle.load(open('fo2d_m6.pkl','rb'))
B = sp.symbols('B6_0:8')
core = []
for c in conds:
    fl = sp.factor_list(c)[1]
    core.append(sp.expand([f for f, mlt in fl if f.has(B[0])][0]))
for i, c in enumerate(core):
    order = max(k for k in range(8) if c.has(B[k]))
    print(f"cond {i}: order {order} in B6; degree in r {sp.Poly(c, r).degree()}")
# express each as sum_k coeff_k(r) B_k
def coeffs(c):
    order = max(k for k in range(8) if c.has(B[k]))
    return [sp.factor(sp.expand(c).coeff(B[k])) for k in range(order+1)]
for i, c in enumerate(core):
    print(f"cond {i} coefficients:", [str(x)[:90] for x in coeffs(c)])
# lowest order condition: solve as ODE
low = min(core, key=lambda c: max(k for k in range(8) if c.has(B[k])))
f = sp.Function('f')(r)
ode = sum(sp.expand(low).coeff(B[k])*sp.diff(f, r, k) for k in range(8) if low.has(B[k]))
print("lowest-order ODE:", sp.factor(ode))
try:
    sol = sp.dsolve(sp.Eq(ode, 0), f)
    print("general solution:", sol)
except Exception as e:
    print("dsolve failed:", e)
pickle.dump(core, open('fo2d_m6_core.pkl','wb'))
