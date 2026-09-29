import sympy as sp, numpy as np, pickle, sys
x, y, z, r, eps = sp.symbols('x y z r epsilon', real=True)
R = sp.sqrt(x**2 + y**2 + z**2); ct = z/R
P = lambda l: sp.legendre(l, ct)

def axis_data(A, B, g_coefs, epsval):
    """Phi = A(r) + eps*B(r)*sum c_l P_l ; return U(z), W(z)=Phi_xx(0,z), Y(z)=Phi_xxxx(0,z) as rational functions of z"""
    g = sum(c*P(l) for l, c in g_coefs.items())
    Phi = A.subs(r, R) + epsval*B.subs(r, R)*g
    U = sp.simplify(Phi.subs({x: 0, y: 0}))
    W = sp.simplify(sp.diff(Phi, x, 2).subs({x: 0, y: 0}))
    Y = sp.simplify(sp.diff(Phi, x, 4).subs({x: 0, y: 0}))
    return [sp.cancel(sp.simplify(e.subs(sp.sqrt(z**2), z))) for e in (U, W, Y)]

if __name__ == '__main__':
    A = -1/(1+r); B = r**2/(1+r)**4
    blind = {2: -sp.Rational(11, 2), 4: -sp.Rational(9, 20), 6: 1}
    for label, gc in [('blind P6-0.45P4-5.5P2', blind), ('pure P2', {2: 1})]:
        U, W, Y = axis_data(A, B, gc, eps)
        print(f"--- {label}")
        print("  W - U'/z =", sp.simplify(W - sp.diff(U, z)/z))
        print("  Y - (spherical Y for U) =", sp.factor(sp.simplify(Y - (3*sp.diff(U, z, 2)/z**2 - 3*sp.diff(U, z)/z**3))))
