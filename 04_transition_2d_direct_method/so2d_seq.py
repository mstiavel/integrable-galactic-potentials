import sympy as sp, pickle, time
r, th, pr, pt = sp.symbols('r theta p_r p_theta', positive=True)
alpha, beta, delta, b = sp.symbols('alpha beta delta b')
A = -1/(1+r**2); Ap = sp.diff(A, r)
B2 = r**2/(1+r**2)**2
H0 = pr**2/2 + pt**2/(2*r**2) + A; I0 = (alpha + beta*H0)*pt**2 + delta*pt**4
import sys
m = int(sys.argv[1]); cm = {2: sp.Rational(3, 8), 0: sp.Rational(1, 4)}[m]
H1r = b*sp.diff(B2, r)*cm; H1th = sp.I*m*b*B2*cm
S = sp.expand(-(-sp.diff(I0, pr)*H1r - sp.diff(I0, pt)*H1th))     # equation: {Q e^{im th}, H0} = S
t0 = time.time()
# unknown coefficient functions
c = {(j, k): sp.Function(f'c{j}{k}')(r) for j in range(5) for k in range(5-j) if (j+k) % 2 == 0}
Q = sum(c[key]*pr**key[0]*pt**key[1] for key in c)
lhs = pr*sp.diff(Q, r) + (pt/r**2)*(sp.I*m)*Q - (Ap - pt**2/r**3)*sp.diff(Q, pr)
eq = sp.expand(lhs - S)
P = sp.Poly(eq, pr, pt)
E = {mon: co for mon, co in zip(P.monoms(), P.coeffs())}
# solve sequentially, top degree first: monomials of degree 5 involve only quartic c's
sol = {}
def solve_ode_first_order(expr, fn):
    """expr linear in fn and fn': a fn' + b fn + s = 0 -> fn = solution with integration constant"""
    expr = sp.expand(expr)
    a_ = expr.coeff(sp.Derivative(fn, r)); rest = sp.expand(expr - a_*sp.Derivative(fn, r))
    b_ = rest.coeff(fn); s_ = sp.expand(rest - b_*fn)
    if a_ == 0:
        return sp.solve(expr, fn)[0]
    mu = sp.exp(sp.integrate(sp.simplify(b_/a_), r))
    mu = sp.simplify(mu)
    K = sp.Symbol(f'K_{fn.func.__name__}')
    integ = sp.integrate(sp.simplify(-s_*mu/a_), r)
    return sp.simplify((integ + K)/mu)
order = [(5, (5, 0), (4, 0)), (5, (4, 1), (3, 1)), (5, (3, 2), (2, 2)), (5, (2, 3), (1, 3)), (5, (1, 4), (0, 4))]
for deg, mon, key in order:
    ex = E[mon].subs(sol).doit()
    sol[c[key]] = solve_ode_first_order(ex, c[key])
    print(f"  c{key} =", sp.simplify(sol[c[key]]), flush=True)
# remaining degree-5 equation: compatibility
comp5 = sp.simplify(E[(0, 5)].subs(sol).doit())
print("  degree-5 compatibility:", comp5, flush=True)
Kc = {n: sp.Symbol(f'K_c{n}') for n in ('40','31','22','13','04','20','11','02','00')}
if m == 2:
    fix = {Kc['31']: 0, Kc['13']: 0, Kc['04']: 0, Kc['40']: 3*b*beta/16}
else:
    fix = {}
sol = {k: sp.simplify(v.subs(fix)) for k, v in sol.items()}
E = {k: v.subs(fix) for k, v in E.items()}
print("  after constraints:", {str(k): v for k, v in sol.items()}, flush=True)
order3 = [((3, 0), (2, 0)), ((2, 1), (1, 1)), ((1, 2), (0, 2))]
for mon, key in order3:
    ex = E[mon].subs(sol).doit()
    sol[c[key]] = solve_ode_first_order(ex, c[key])
    print(f"  c{key} =", str(sp.simplify(sol[c[key]]))[:400], flush=True)
comp3 = sp.simplify(E[(0, 3)].subs(sol).doit())
print("  degree-3 compatibility:", str(comp3)[:600], flush=True)
ex = E[(1, 0)].subs(sol).doit(); sol[c[(0, 0)]] = solve_ode_first_order(ex, c[(0, 0)])
print("  c00 =", str(sp.simplify(sol[c[(0, 0)]]))[:400], flush=True)
comp1 = sp.simplify(E[(0, 1)].subs(sol).doit())
print("  degree-1 compatibility:", str(comp1)[:600], flush=True)
print("time", time.time()-t0)
pickle.dump({str(k): v for k, v in sol.items()}, open(f'cored_mode{m}_seq.pkl', 'wb'))
