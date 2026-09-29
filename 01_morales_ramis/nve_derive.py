import sympy as sp

x, y, z, r, E = sp.symbols('x y z r E', positive=True)
A = sp.Function('A')
B = sp.Function('B')

R = sp.sqrt(x**2 + y**2 + z**2)
P2 = (3*z**2/R**2 - 1)/2
Phi = A(R) + B(R)*P2

def hess(i, j):
    return sp.diff(Phi, i, j)

# --- z-axis orbit: point (0,0,z). transverse second derivative d2Phi/dx2
Wax = sp.simplify(hess(x, x).subs({x: 0, y: 0}).subs(z, r))
Wax_claim = (sp.diff(A(r), r) + sp.diff(B(r), r))/r - 3*B(r)/r**2
print("axis W - claim:", sp.simplify(Wax - Wax_claim))
# potential along axis
Uax = sp.simplify(Phi.subs({x: 0, y: 0}).subs(z, r))
print("axis U:", Uax)

# --- x-axis (equatorial rectilinear, Lz=0) orbit: point (x,0,0). transverse d2Phi/dz2 (out of plane) and d2Phi/dy2 (in plane)
Weq_z = sp.simplify(hess(z, z).subs({y: 0, z: 0}).subs(x, r))
Weq_y = sp.simplify(hess(y, y).subs({y: 0, z: 0}).subs(x, r))
Veq = sp.simplify(Phi.subs({y: 0, z: 0}).subs(x, r))
Weq_z_claim = (sp.diff(A(r), r) - sp.diff(B(r), r)/2)/r + 3*B(r)/r**2
print("eq  W_z - claim:", sp.simplify(Weq_z - Weq_z_claim))
print("eq  W_y:", Weq_y, "  (= V'/r, the central-force in-plane NVE)")
print("eq  V:", Veq)

# --- check the NVE in z-variable form: 2(E-U) xi'' - U' xi' + W xi = 0
# with spherical case B=0: xi=z must be exact solution
U = A(r); W = sp.diff(U, r)/r
xi = r
print("spherical: residual of xi=r:", sp.simplify(2*(E-U)*sp.diff(xi, r, 2) - sp.diff(U, r)*sp.diff(xi, r) + W*xi))

# --- indicial exponents at a power-law singular point: U ~ u r^mu, B ~ b r^mu
u, b, mu, rho = sp.symbols('u b mu rho')
Ul = u*r**mu; Bl = b*r**mu
Wl = sp.diff(Ul, r)/r - 3*Bl/r**2
expr = sp.expand(sp.simplify((2*(-Ul)*sp.diff(r**rho, r, 2) - sp.diff(Ul, r)*sp.diff(r**rho, r) + Wl*r**rho)/r**(rho+mu-2)))
print("indicial polynomial:", sp.factor(expr))
print("roots:", sp.solve(expr, rho))
