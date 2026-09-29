import sympy as sp, pickle
from sf_general import Fj, S, C, redSC
M, rhs, u0, pool, usol, conds, params = pickle.load(open('deg6_d4_mu1.pkl','rb'))
a0, a1 = params; K = sp.symbols('K')
g = a0*S**2 + a1*C**2
th = sp.symbols('theta')
# cond1 => f'' = K sqrt(g) - f ; cond3 => g(389 f - 241 f'') = 630 (a0-a1) S C f'
c3 = [fa for fa, m in sp.factor_list(conds[2])[1] if fa.has(Fj[0])][0]
c2 = [fa for fa, m in sp.factor_list(conds[1])[1] if fa.has(Fj[0])][0]
sq = sp.sqrt(g)
F2e = K*sq - Fj[0]
F1e = sp.solve(c3.subs(Fj[2], F2e), Fj[1])[0]
def D(e): return sp.diff(e, S)*C - sp.diff(e, C)*S + sp.diff(e, Fj[0])*Fj[1]
# f'' from differentiating F1e, with f' replaced, must equal F2e
alg = sp.simplify(D(F1e).subs(Fj[1], F1e) - F2e)
print("algebraic equation for f (numerator):", sp.factor(sp.numer(sp.together(alg))))
fsol = sp.solve(alg, Fj[0])
print("f =", fsol)
# check consistency: f' from F1e must equal derivative of fsol
for fs in fsol:
    fp = sp.simplify(sp.diff(fs, S)*C - sp.diff(fs, C)*S)
    resid = sp.simplify(fp - F1e.subs(Fj[0], fs))
    print("  residual of cond3-first-order on this f:", sp.factor(sp.numer(sp.together(resid))))
    # cond 2 (fourth order): substitute derivatives
    derivs = {Fj[0]: fs}
    cur = fs
    for k in range(1, 5):
        cur = sp.simplify(sp.diff(cur, S)*C - sp.diff(cur, C)*S); derivs[Fj[k]] = cur
    r2 = sp.simplify(c2.subs(derivs))
    print("  cond2 residual:", sp.factor(sp.numer(sp.together(r2))))
