import sympy as sp, sys, time
exec(open('sf_w6.py').read().split("print(\"\\n--- Q relations directly")[0].replace("w = int(sys.argv[1])", "w = 2"))
a, b, c0 = params
rel = sp.factor(redSC(al['b20']*be['b11'] - al['b11']*be['b20']))
# pick the factor that is linear in (F0,F1,F2) and contains the parameters
fac = [fa for fa, m in sp.factor_list(rel)[1] if fa.has(Fj[1]) and fa.has(a) ]
master = fac[0]
print("master factor:", master)
F2sol = sp.solve(master, Fj[2])[0]
def D2(e): return sp.diff(e, S)*C - sp.diff(e, C)*S + sp.diff(e, Fj[0])*Fj[1] + sp.diff(e, Fj[1])*Fj[2]
F3sol = sp.cancel(D2(F2sol).subs(Fj[2], F2sol))
print("\nresiduals of the Q's on the master ODE:")
for n, q in Q.items():
    v = sp.cancel(sp.together(q.subs(Fj[3], F3sol).subs(Fj[2], F2sol)))
    v = sp.factor(redSC(sp.numer(v)))
    print(" ", n, ":", str(v)[:400])
print("\n--- solving for Killing parameters that kill the residuals identically in theta:")
conds = []
for n, q in Q.items():
    v = sp.cancel(sp.together(q.subs(Fj[3], F3sol).subs(Fj[2], F2sol)))
    v = redSC(sp.numer(v))
    # remove the F1 and trigonometric prefactors: collect coefficients in S, C and F's
    P = sp.Poly(v, S, C, *Fj[:2])
    conds += [co for co in P.coeffs()]
sol = sp.solve(conds, [a, b, c0], dict=True)
print("solutions:", sol)
