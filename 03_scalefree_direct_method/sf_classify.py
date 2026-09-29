import sympy as sp, sys, pickle, time
from sf_direct import system, th, f
J = 8
Fj = sp.symbols(f'F0:{J}')          # f and derivatives
Bj = {n: sp.symbols(f'{n}_0:{J}') for n in ('b20', 'b11', 'b02', 'c')}

def jets(eqs, fns):
    b20, b11, b02, c = fns
    rep = {}
    for fn, name in zip(fns, ('b20', 'b11', 'b02', 'c')):
        for k in range(J-1, -1, -1):
            rep[sp.Derivative(fn, (th, k)) if k else fn] = Bj[name][k]
    for k in range(J-1, -1, -1):
        rep[sp.Derivative(f, (th, k)) if k else f] = Fj[k]
    return {key: sp.expand(e.doit().subs(rep)) for key, e in eqs}

def Dth(e):
    out = sp.diff(e, th)
    for name in Bj:
        for k in range(J-1): out += sp.diff(e, Bj[name][k])*Bj[name][k+1]
    for k in range(J-1): out += sp.diff(e, Fj[k])*Fj[k+1]
    return sp.expand(out)

def classify(w=4, muval=1, params=None):
    t0 = time.time()
    if params is None: params = sp.symbols('a b')
    eqs, fns = system(w, params, muval)
    E = jets(eqs, fns)
    E1, E2, E3, E4, E5, E6 = E[(3,3,0)], E[(3,2,1)], E[(3,1,2)], E[(3,0,3)], E[(1,1,0)], E[(1,0,1)]
    b20p = sp.solve(E1, Bj['b20'][1])[0]
    b02p = sp.solve(E4, Bj['b02'][1])[0]
    b11p = sp.solve(E2.subs(Bj['b20'][1], b20p), Bj['b11'][1])[0]
    subsP = {Bj['b20'][1]: b20p, Bj['b02'][1]: b02p, Bj['b11'][1]: b11p}
    def red(e):   # replace first derivatives of b's by their expressions, repeatedly (for second derivatives use Dth then substitute)
        return sp.expand(e.subs(subsP))
    R1 = sp.numer(sp.together(red(E3)))
    c0, c1 = Bj['c'][0], Bj['c'][1]
    T5 = -sp.expand(E5.subs({c0: 0, c1: 0})); T6 = -sp.expand(E6.subs({c0: 0, c1: 0}))
    # check structure: E5 = cos c' + 4 sin c - T5 ; E6 = -sin c' + 4 cos c - T6  (w=4, mu=1: weight factor 4)
    csol = sp.expand((sp.sin(th)*T5 + sp.cos(th)*T6)/w)
    cp = sp.expand(sp.cos(th)*T5 - sp.sin(th)*T6)
    assert sp.simplify(E5 - (sp.cos(th)*c1 + w*sp.sin(th)*c0 - T5)) == 0
    assert sp.simplify(E6 - (-sp.sin(th)*c1 + w*sp.cos(th)*c0 - T6)) == 0
    R3 = sp.numer(sp.together(red(Dth(csol) - cp)))
    R3 = sp.expand(red(R3))
    R1p = sp.expand(red(Dth(R1)))
    print(f"  relations built in {time.time()-t0:.0f}s; sizes {len(str(R1))}, {len(str(R3))}, {len(str(R1p))}", flush=True)
    Bvars = [Bj['b20'][0], Bj['b11'][0], Bj['b02'][0]]
    # all three relations are linear in Bvars
    M, rhs = sp.linear_eq_to_matrix([R1, R3, R1p], Bvars)
    det = sp.factor(M.det())
    print("  det of 3x3 system for (b20,b11,b02):", det, flush=True)
    return dict(M=M, rhs=rhs, R1=R1, R3=R3, R1p=R1p, red=red, subsP=subsP, csol=csol, cp=cp, Bvars=Bvars, params=params)

if __name__ == '__main__':
    res = classify()
    pickle.dump({k: v for k, v in res.items() if k not in ('red',)}, open('classify_w4.pkl', 'wb'))
