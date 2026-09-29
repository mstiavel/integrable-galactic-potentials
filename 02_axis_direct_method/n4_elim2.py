import sympy as sp, pickle, sys
z = sp.symbols('z'); alpha, beta = sp.symbols('alpha beta')
out = pickle.load(open('n4_first_order.pkl','rb'))
B1 = sp.Function('B1')(z)
G = sp.symbols('G0:5'); Bs = sp.symbols('Bd0:8')

def Dz(e):
    return sp.diff(e, z) + sum(sp.diff(e, G[k])*G[k+1] for k in range(4)) + sum(sp.diff(e, Bs[k])*Bs[k+1] for k in range(7))

def get_ode(orbit, av, bv):
    g, ks, e3, e1 = out[orbit]
    def tosym(e):
        e = e.subs({alpha: av, beta: bv}); rep = {}
        for k in range(4, -1, -1): rep[sp.Derivative(g, (z, k)) if k else g] = G[k]
        for k in range(7, -1, -1): rep[sp.Derivative(B1, (z, k)) if k else B1] = Bs[k]
        return sp.together(e.subs(rep))
    E3 = sp.numer(tosym(e3)); E1 = sp.numer(tosym(e1))
    G3_E3 = sp.solve(E3, G[3])[0]
    G2_E1 = sp.solve(E1, G[2])[0]
    R1 = sp.numer(sp.together(Dz(G2_E1).subs(G[2], G2_E1) - G3_E3))
    R1 = sp.expand(R1)
    a1, a0 = R1.coeff(G[1]), R1.coeff(G[0])
    S = sp.expand(R1 - a1*G[1] - a0*G[0])
    assert not S.has(G[0]) and not S.has(G[1])
    p = sp.cancel(-a0/a1); s = sp.cancel(-S/a1)
    # consistency of gamma'' from E1 with derivative of gamma' = p gamma + s
    L1 = G2_E1.subs(G[1], p*G[0] + s)
    cons = sp.expand(sp.numer(sp.together(L1 - ((sp.diff(p, z) + p**2)*G[0] + p*s + Dz(s)))))
    cg = cons.coeff(G[0])
    T = sp.expand(cons - cg*G[0])
    return sp.factor(cg), T, p, s, ks

if __name__ == '__main__':
    av, bv = sp.Rational(1), sp.Rational(1, 3)
    res = {}
    for orbit in ['axis','eq']:
        cg, T, p, s, ks = get_ode(orbit, av, bv)
        print(f"\n{orbit}: coefficient of gamma in consistency:", cg)
        order = max(k for k in range(8) if T.has(Bs[k]))
        print(f"{orbit}: ODE for B1: order {order}, inhomogeneous in", [k for k in ks if T.has(k)])
        # leading structure: coefficient polynomials of B derivatives
        for k in range(order, -1, -1):
            print(f"   coeff of B^({k}):", sp.factor(T.coeff(Bs[k])))
        print("   inhomogeneous part:", sp.factor(T.subs({b: 0 for b in Bs})))
        res[orbit] = (T, p, s, ks)
        sys.stdout.flush()
    pickle.dump(res, open('n4_odes.pkl','wb'))
