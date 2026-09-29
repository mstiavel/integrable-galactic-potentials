import numpy as np
from scipy.integrate import solve_ivp
import itertools, sys

# Family: Phi = A(r) + B(r) P2(cos theta)
#   A = -1/(1+r)            Hernquist cusp (rho ~ r^-1), spherical centre
#   B = eps r^2/(1+r)^n     quadrupole, ~ r^2 at centre (spherical centre), decays outward
# Axis orbit NVE (independent variable z):  2(E-U) xi'' - U' xi' + W xi = 0
#   U = A+B,  W = U'/z - 3B/z^2
# Equatorial orbit NVE (out of plane):     2(E-V) xi'' - V' xi' + W xi = 0
#   V = A-B/2, W = V'/x + 3B/x^2

def make(eps, n, orbit):
    def A(z):  return -1/(1+z)
    def dA(z): return 1/(1+z)**2
    def Bf(z): return eps*z**2/(1+z)**n
    def dB(z): return eps*(2*z/(1+z)**n - n*z**2/(1+z)**(n+1))
    if orbit == 'axis':
        U  = lambda z: A(z) + Bf(z)
        dU = lambda z: dA(z) + dB(z)
        W  = lambda z: dU(z)/z - 3*Bf(z)/z**2
    else:
        U  = lambda z: A(z) - Bf(z)/2
        dU = lambda z: dA(z) - dB(z)/2
        W  = lambda z: dU(z)/z + 3*Bf(z)/z**2
    return U, dU, W

def turning_points(eps, n, E, orbit):
    # numerator of E-U as polynomial in z, computed numerically via sampling/fit is fragile;
    # do it exactly: E - U = [E(1+z)^n + s*(1+z)^(n-1) - c*eps z^2]/(1+z)^n  with (s,c)=(1,1) axis, (1,-1/2) eq
    s = 1.0; c = 1.0 if orbit == 'axis' else -0.5
    p = E*np.poly1d([1, 1])**n + s*np.poly1d([1, 1])**(n-1) - c*eps*np.poly1d([1, 0, 0])
    return np.roots(p.coeffs)

def integrate_path(U, dU, W, E, pts, rtol=1e-11, atol=1e-13):
    """Propagate the 2x2 fundamental matrix (basis xi, xi') along a polygonal path of complex points."""
    M = np.eye(2, dtype=complex)
    for z0, z1 in zip(pts[:-1], pts[1:]):
        dz = z1 - z0
        def rhs(s, y):
            zz = z0 + s*dz
            Y = y.reshape(2, 2)
            # d/dz [xi, xi'] = [xi', (U' xi' - W xi)/(2(E-U))]
            a = dU(zz)/(2*(E-U(zz))); b = -W(zz)/(2*(E-U(zz)))
            D = np.array([[0, 1], [b, a]], dtype=complex)
            return (D @ Y * dz).reshape(-1)
        sol = solve_ivp(rhs, (0, 1), M.reshape(-1), method='DOP853', rtol=rtol, atol=atol)
        M = sol.y[:, -1].reshape(2, 2)
    return M

def loop(center, base, sing, nseg=64):
    """Polygonal loop from base around 'center' only: radial approach, circle, return."""
    others = [s for s in sing if abs(s-center) > 1e-9]
    rad = 0.4*min([abs(center-s) for s in others] + [abs(center-base)])
    direction = (base-center)/abs(base-center)
    start = center + rad*direction
    # approach: straight from base to start, but detour if passing near another singular point
    approach = [base, start]
    th0 = np.angle(direction)
    circle = [center + rad*np.exp(1j*(th0 + 2*np.pi*k/nseg)) for k in range(nseg+1)]
    return approach + circle[1:] + [start, base]

def sl2(M):
    d = np.linalg.det(M)
    return M/np.sqrt(d)

def comm(P, Q):
    return P @ Q @ np.linalg.inv(P) @ np.linalg.inv(Q)

def run(eps, n, E, orbit, base=0.6+0.9j):
    U, dU, W = make(eps, n, orbit)
    tps = turning_points(eps, n, E, orbit)
    sing = [0.0+0j, -1.0+0j] + list(tps)
    print(f"\n=== orbit={orbit}  eps={eps}  n={n}  E={E}")
    print("singular points:", np.round(sing, 4))
    Ms = []
    for s in sing:
        pts = loop(s, base, sing)
        M = integrate_path(U, dU, W, E, pts)
        Ms.append(M)
        ev = np.linalg.eigvals(M)
        print(f"  around {s:.4f}: det={np.linalg.det(M):.6f}  eig={np.round(ev,5)}")
    # consistency: product of all loops (suitably ordered) ~ loop around infinity; skip.
    # derived subgroup test
    S = [sl2(M) for M in Ms]
    comms = []
    for i, j in itertools.combinations(range(len(S)), 2):
        C = comm(S[i], S[j])
        comms.append(C)
    # commutators mutually commuting?  (necessary for virtually-abelian G)
    worst = 0.0
    for C1, C2 in itertools.combinations(comms, 2):
        worst = max(worst, np.linalg.norm(comm(C1, C2) - np.eye(2)))
    print(f"  max ||[C,C']-I|| over commutator pairs : {worst:.3e}")
    # infinite-order element? look at traces of products
    hyper = 0.0
    for i, j in itertools.permutations(range(len(S)), 2):
        t = abs(np.trace(S[i] @ S[j]))
        hyper = max(hyper, t)
    for i, j, k in itertools.permutations(range(len(S)), 3):
        t = abs(np.trace(S[i] @ S[j] @ S[k]))
        hyper = max(hyper, t)
    print(f"  max |tr| over short products (in SL2 normalisation): {hyper:.4f}  (>2 => hyperbolic, infinite order)")
    return worst, hyper

if __name__ == '__main__':
    n = 4
    for orbit in ['axis', 'eq']:
        run(0.0, n, -0.5, orbit)      # spherical control
        run(0.3, n, -0.5, orbit)
        run(1.0, n, -0.5, orbit)
        run(1.0, n, -0.2, orbit)
