# First-order n=4 direct-method test with two multipoles (P2, P4) around a spherical Hernquist cusp.
# Verifies: axis condition = P2-only operator on (3 B2 + 10 B4)/3, equator = -(same) on (3 B2 - 7.5 B4)/3.
# Requires direct_axis.py and the pickle n4_odes.pkl produced by n4_elim2.py in the same directory.
import sympy as sp, pickle
from direct_axis import hierarchy
zz, eps = sp.symbols('z epsilon'); alpha, beta = sp.symbols('alpha beta')
Uf = sp.Function('U')(zz); Wf = sp.Function('W')(zz)
rem, cs, sol = hierarchy(4, 0); c0f, c2f = cs[0], cs[2]; E3 = rem[1][1]; E1 = rem[2][1]
b2 = sp.Function('b2')(zz); b4 = sp.Function('b4')(zz)
G = sp.symbols('G0:5'); Bs2 = sp.symbols('P0:8'); Bs4 = sp.symbols('Q0:8')
def Dz(e): return sp.diff(e, zz) + sum(sp.diff(e, G[k])*G[k+1] for k in range(4)) + sum(sp.diff(e, Bs2[k])*Bs2[k+1] for k in range(7)) + sum(sp.diff(e, Bs4[k])*Bs4[k+1] for k in range(7))
A0 = -1/(1+zz)
def ode_for(orbit, av, bv):
    g = sp.Function('g')(zz); k0, k1, k2 = sp.symbols('k0 k1 k2')
    if orbit == 'axis':
        U = A0 + eps*(b2 + b4); W = sp.diff(U, zz)/zz - eps*(3*b2 + 10*b4)/zz**2
    else:
        U = A0 + eps*(-b2/2 + sp.Rational(3,8)*b4); W = sp.diff(U, zz)/zz + eps*(3*b2 - sp.Rational(15,2)*b4)/zz**2
    c0 = (alpha + beta*A0)*zz**2 + eps*g; c2 = beta*zz**2/2 + eps*(k0 + k1*zz + k2*zz**2)
    subs = {c0f: c0, c2f: c2, Uf: U, Wf: W}
    e3 = sp.diff(sp.expand(E3.subs(subs).doit()), eps).subs(eps, 0); e1 = sp.diff(sp.expand(E1.subs(subs).doit()), eps).subs(eps, 0)
    def tosym(e):
        e = e.subs({alpha: av, beta: bv}); rep = {}
        for k in range(4, -1, -1): rep[sp.Derivative(g, (zz, k)) if k else g] = G[k]
        for k in range(7, -1, -1):
            rep[sp.Derivative(b2, (zz, k)) if k else b2] = Bs2[k]; rep[sp.Derivative(b4, (zz, k)) if k else b4] = Bs4[k]
        return sp.together(e.subs(rep))
    E3s = sp.numer(tosym(e3)); E1s = sp.numer(tosym(e1))
    G3 = sp.solve(E3s, G[3])[0]; G2 = sp.solve(E1s, G[2])[0]
    R1 = sp.expand(sp.numer(sp.together(Dz(G2).subs(G[2], G2) - G3)))
    a1, a0 = R1.coeff(G[1]), R1.coeff(G[0]); S = sp.expand(R1 - a1*G[1] - a0*G[0])
    p = sp.cancel(-a0/a1); s = sp.cancel(-S/a1)
    L1 = G2.subs(G[1], p*G[0] + s)
    cons = sp.expand(sp.numer(sp.together(L1 - ((sp.diff(p, zz) + p**2)*G[0] + p*s + Dz(s)))))
    assert cons.coeff(G[0]) == 0
    return sp.expand(cons - cons.coeff(G[0])*G[0])
if __name__ == '__main__':
    av, bv = sp.Rational(1), sp.Rational(1,3)
    Tax = ode_for('axis', av, bv); Teq = ode_for('eq', av, bv)
    zero_k = {sp.Symbol('k0'):0, sp.Symbol('k1'):0, sp.Symbol('k2'):0}
    print("axis ODE vanishes on 3b2+10b4=0:", sp.simplify(Tax.subs(zero_k).subs({Bs4[k]: -sp.Rational(3,10)*Bs2[k] for k in range(8)})) == 0)
    print("equator ODE vanishes on 3b2-7.5b4=0:", sp.simplify(Teq.subs(zero_k).subs({Bs4[k]: sp.Rational(2,5)*Bs2[k] for k in range(8)})) == 0)
