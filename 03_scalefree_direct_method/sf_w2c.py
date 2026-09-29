import sympy as sp, sys, time, pickle
exec(open('sf_w6.py').read().split("print(\"\\n--- Q relations directly")[0].replace("w = int(sys.argv[1])", "w = 2"))
a, b, c0 = params
rel = sp.factor(redSC(al['b20']*be['b11'] - al['b11']*be['b20']))
fac = [fa for fa, m in sp.factor_list(rel)[1] if fa.has(Fj[1]) and fa.has(a)]
master = fac[0]
F2sol = sp.solve(master, Fj[2])[0]
def D2(e): return sp.diff(e, S)*C - sp.diff(e, C)*S + sp.diff(e, Fj[0])*Fj[1] + sp.diff(e, Fj[1])*Fj[2]
F3sol = sp.cancel(D2(F2sol).subs(Fj[2], F2sol))
conds = set()
for n, q in Q.items():
    v = sp.cancel(sp.together(q.subs(Fj[3], F3sol).subs(Fj[2], F2sol)))
    v = redSC(sp.numer(v))
    P = sp.Poly(v, S, C, Fj[0], Fj[1])
    for co in P.coeffs():
        co = sp.factor(co)
        conds.add(co)
conds = list(conds)
print("number of coefficient conditions:", len(conds), flush=True)
G = sp.groebner(conds, a, b, c0, order='lex')
print("Groebner basis:", G, flush=True)
print("solutions:", sp.solve(list(G), [a, b, c0], dict=True))
pickle.dump((master, list(G)), open('w2_result.pkl','wb'))
