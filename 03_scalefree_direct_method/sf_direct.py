import sympy as sp, itertools, sys
th, r = sp.symbols('theta r', positive=True)
px, pz = sp.symbols('p_x p_z')
x = r*sp.sin(th); z = r*sp.cos(th)
f = sp.Function('f')(th)
mu = sp.Symbol('mu')
V = r**mu*f
L = x*pz - z*px

def quartic_A(w, params):
    if w == 2: a, b, c = params; return a*px**4 + b*px**2*pz**2 + c*pz**4
    if w == 4: a, b = params;   return L**2*(a*px**2 + b*pz**2)
    if w == 6: (a,) = params;   return a*L**4

def dq(expr, which):
    # partial derivative w.r.t. Cartesian x or z of an expression in (r, theta)
    if which == 'x': return sp.sin(th)*sp.diff(expr, r) + sp.cos(th)/r*sp.diff(expr, th)
    else:            return sp.cos(th)*sp.diff(expr, r) - sp.sin(th)/r*sp.diff(expr, th)

def system(w, params, muval=1):
    A = quartic_A(w, params)
    b20, b11, b02, c = [sp.Function(n)(th) for n in ('b20', 'b11', 'b02', 'c')]
    d = w - 2
    B = r**(d+muval)*(b20*px**2 + b11*px*pz + b02*pz**2)
    C = r**(d+2*muval)*c
    Vm = V.subs(mu, muval)
    Vx, Vz = dq(Vm, 'x'), dq(Vm, 'z')
    # cubic part:  p.grad_q B - grad_p A . grad V = 0
    cub = px*dq(B, 'x') + pz*dq(B, 'z') - (sp.diff(A, px)*Vx + sp.diff(A, pz)*Vz)
    # linear part: p.grad_q C - grad_p B . grad V = 0
    lin = px*dq(C, 'x') + pz*dq(C, 'z') - (sp.diff(B, px)*Vx + sp.diff(B, pz)*Vz)
    eqs = []
    for expr, deg in ((cub, 3), (lin, 1)):
        P = sp.Poly(sp.expand(expr), px, pz)
        for (i, j), co in zip(P.monoms(), P.coeffs()):
            co = sp.simplify(co / r**(d+muval-1 if deg == 3 else d+2*muval-1))
            co = sp.powsimp(sp.expand(co), force=True)
            assert not co.has(r), (deg, i, j, co)
            eqs.append(((deg, i, j), sp.simplify(co)))
    return eqs, (b20, b11, b02, c)

if __name__ == '__main__':
    for w, params in [(4, sp.symbols('a b')), (2, sp.symbols('a b c')), (6, sp.symbols('a',))]:
        params = tuple(params) if isinstance(params, (tuple, list)) else (params,)
        eqs, fns = system(w, params)
        print(f"\n===== weight w={w}, Killing parameters {params}")
        for key, e in eqs:
            print(f"  [deg {key[0]}, px^{key[1]} pz^{key[2]}]:  {e} = 0")
