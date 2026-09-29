import sympy as sp, numpy as np, itertools, sys
from scipy.integrate import solve_ivp
from numpy.polynomial import polynomial as Pn
from ve3_setup import axis_data
z, r, eps = sp.symbols('z r epsilon', real=True)

def build_system(U, W, Y, E):
    Uz, Wz, Yz = [sp.lambdify(z, e, 'numpy') for e in (U, W, Y)]
    dU = sp.lambdify(z, sp.diff(U, z), 'numpy'); ddU = sp.lambdify(z, sp.diff(U, z, 2), 'numpy')
    dW = sp.lambdify(z, sp.diff(W, z), 'numpy')
    def coef(zz):
        den = 2*(E - Uz(zz))
        p1 = -Wz(zz)/den; q1 = dU(zz)/den; p2 = -ddU(zz)/den
        s2 = -dW(zz)/(2*den); s3 = -dW(zz)/den; t3 = -Yz(zz)/(6*den)
        return p1, q1, p2, s2, s3, t3
    # monomial index: 0 xi,1 xi',2 xi^2,3 xi xi',4 xi'^2,5 ze,6 ze',7 xi ze,8 xi ze',9 xi' ze,10 xi' ze',11 xi^3,12 xi^2xi',13 xi xi'^2,14 xi'^3,15 xi3,16 xi3'
    def Amat(zz):
        p1, q1, p2, s2, s3, t3 = coef(zz)
        A = np.zeros((17, 17), dtype=complex)
        A[0, 1] = 1; A[1, 0] = p1; A[1, 1] = q1
        A[2, 3] = 2; A[3, 2] = p1; A[3, 3] = q1; A[3, 4] = 1; A[4, 3] = 2*p1; A[4, 4] = 2*q1
        A[5, 6] = 1; A[6, 5] = p2; A[6, 6] = q1; A[6, 2] = s2
        A[7, 9] = 1; A[7, 8] = 1
        A[8, 7] = p2; A[8, 8] = q1; A[8, 10] = 1; A[8, 11] = s2
        A[9, 7] = p1; A[9, 9] = q1; A[9, 10] = 1
        A[10, 8] = p1; A[10, 9] = p2; A[10, 10] = 2*q1; A[10, 12] = s2
        A[11, 12] = 3; A[12, 11] = p1; A[12, 12] = q1; A[12, 13] = 2
        A[13, 12] = 2*p1; A[13, 13] = 2*q1; A[13, 14] = 1; A[14, 13] = 3*p1; A[14, 14] = 3*q1
        A[15, 16] = 1; A[16, 15] = p1; A[16, 16] = q1; A[16, 7] = s3; A[16, 11] = t3
        return A
    return Amat

def singular_points(U, W, Y, E):
    num, den = sp.fraction(sp.cancel(E - U))
    pts = [complex(p) for p in sp.Poly(num, z).nroots(n=30)]
    for e in (U, W, Y):
        d = sp.fraction(sp.cancel(e))[1]
        pd = sp.Poly(sp.sqf_part(d), z)
        if pd.degree() > 0: pts += [complex(p) for p in pd.nroots(n=30)]
    out = []
    for p in pts:
        if all(abs(p - q) > 1e-7 for q in out): out.append(p)
    return out

def monodromies(Amat, sing, base=0.6+0.9j, frac=0.4, nseg=128, rtol=1e-11, atol=1e-13):
    Ms = []
    for s in sing:
        others = [t for t in sing if abs(t - s) > 1e-9]
        rad = frac*min([abs(s - t) for t in others] + [abs(s - base)])
        d = (base - s)/abs(base - s); start = s + rad*d; th0 = np.angle(d)
        circ = [s + rad*np.exp(1j*(th0 + 2*np.pi*k/nseg)) for k in range(nseg + 1)]
        pts = [base, start] + circ[1:] + [start, base]
        M = np.eye(17, dtype=complex)
        for z0, z1 in zip(pts[:-1], pts[1:]):
            dz = z1 - z0
            def rhs(t, y):
                Yv = y.reshape(17, 17); return (Amat(z0 + t*dz) @ Yv * dz).reshape(-1)
            sol = solve_ivp(rhs, (0, 1), M.reshape(-1), method='DOP853', rtol=rtol, atol=atol)
            M = sol.y[:, -1].reshape(17, 17)
        Ms.append(M)
    return Ms

def comm_measure(P, Q):
    return np.linalg.norm(P @ Q - Q @ P)/(np.linalg.norm(P)*np.linalg.norm(Q))

def analyse(Ms, sing, label):
    I = np.eye(17)
    # elements of the identity component: unipotent monodromies and squares of finite-order ones, plus conjugates
    elems = []
    for M in Ms:
        ev = np.linalg.eigvals(M)
        if np.allclose(ev, 1, atol=1e-6): elems.append(M)             # unipotent
        else: elems.append(M @ M)                                     # square (eigenvalues +-1 -> unipotent)
    conj = []
    for M in Ms:
        Minv = np.linalg.inv(M)
        for U in elems: conj.append(M @ U @ Minv)
    allel = elems + conj
    # exclude near-identity elements
    allel = [U for U in allel if np.linalg.norm(U - I) > 1e-8]
    worst = 0.0
    for P, Q in itertools.combinations(allel, 2):
        worst = max(worst, comm_measure(P, Q))
    unip = [np.allclose(np.linalg.eigvals(U), 1, atol=1e-5) for U in allel]
    print(f"  {label}: {len(allel)} G0-elements (all unipotent: {all(unip)}), max relative commutator {worst:.2e}")
    return worst

if __name__ == '__main__':
    A = -1/(1+r); B = r**2/(1+r)**4
    blind = {2: -sp.Rational(11, 2), 4: -sp.Rational(9, 20), 6: 1}
    E = -0.5
    for label, gc, ev in [('spherical control (eps=0)', blind, 0), ('blind class, eps=0.05', blind, sp.Rational(1, 20)), ('blind class, eps=0.2', blind, sp.Rational(1, 5)), ('pure P2 control, eps=0.2', {2: 1}, sp.Rational(1, 5))]:
        U, W, Y = axis_data(A, B, gc, ev)
        sing = singular_points(U, W, Y, E)
        Amat = build_system(U, W, Y, E)
        Ms = monodromies(Amat, sing)
        Pinf = np.eye(17, dtype=complex)
        order = np.argsort([np.angle(s - (0.6+0.9j)) for s in sing])
        for k in order: Pinf = Ms[k] @ Pinf
        print(f"--- {label}: singular points {np.round(sing, 3)}; |tr(loop at inf) - 17| = {abs(np.trace(Pinf) - 17):.2e}")
        analyse(Ms, sing, label)
        sys.stdout.flush()
