import sympy as sp, numpy as np, sys, time, itertools
from scipy.signal import fftconvolve

R, z, x, t = sp.symbols('R z x t', real=True)

def taylor_V(Phi, R0, N):
    """Taylor coefficients of V(x,z) = Lz^2/(2(R0+x)^2) + Phi(R0+x, z) about (0,0), total degree <= N, as dict (i,j)->coef"""
    PhiR = sp.diff(Phi, R)
    Lz2 = R0**3*PhiR.subs({R: R0, z: 0})
    V = Lz2/(2*(R0 + t*x)**2) + Phi.subs({R: R0 + t*x, z: t*z})
    ser = sp.series(V, t, 0, N+1).removeO()
    ser = sp.expand(ser)
    coefs = {}
    P = sp.Poly(ser, t, x, z)
    for (k, i, j), c in zip(P.monoms(), P.coeffs()):
        coefs[(i, j)] = complex(sp.N(c, 20))
    return coefs, float(sp.N(Lz2))

class Poly4:
    """dense polynomial in (z1, zb1, z2, zb2), truncated at total degree N"""
    def __init__(self, N):
        self.N = N; self.a = np.zeros((N+1,)*4, dtype=complex)
    def copy(self):
        p = Poly4(self.N); p.a = self.a.copy(); return p
    def truncate(self):
        idx = np.indices(self.a.shape).sum(axis=0)
        self.a[idx > self.N] = 0
        return self
    def degree_slice(self, n):
        idx = np.indices(self.a.shape).sum(axis=0)
        out = Poly4(self.N); out.a = np.where(idx == n, self.a, 0); return out
    def maxabs(self, n):
        idx = np.indices(self.a.shape).sum(axis=0)
        v = np.abs(self.a[idx == n]); return v.max() if v.size else 0.0

def mul(p, q):
    out = Poly4(p.N)
    full = fftconvolve(p.a, q.a)
    out.a = full[:p.N+1, :p.N+1, :p.N+1, :p.N+1]
    return out.truncate()

def deriv(p, axis):
    out = Poly4(p.N)
    n = p.N + 1
    sl = [slice(None)]*4; sl[axis] = slice(1, None)
    idx = np.arange(1, n)
    shape = [1, 1, 1, 1]; shape[axis] = n - 1
    src = p.a[tuple(sl)] * idx.reshape(shape)
    dst = [slice(None)]*4; dst[axis] = slice(0, n-1)
    out.a[tuple(dst)] = src
    return out

def bracket(f, g):
    """{f,g} = -i sum_k (f_z g_zb - f_zb g_z), axes: 0=z1,1=zb1,2=z2,3=zb2"""
    out = Poly4(f.N)
    for zk, zbk in ((0, 1), (2, 3)):
        out.a += mul(deriv(f, zk), deriv(g, zbk)).a - mul(deriv(f, zbk), deriv(g, zk)).a
    out.a *= -1j
    return out.truncate()

def to_complex_vars(coefs, kappa, nu, N):
    """V(x,z) + kinetic -> polynomial in complex variables. x = (z1+zb1)/sqrt(2 kappa), pR = -i sqrt(kappa/2)(z1 - zb1), etc."""
    H = Poly4(N)
    # build x, z, pR, pz as Poly4
    def lin(coefs4):
        p = Poly4(N)
        for axis, c in enumerate(coefs4):
            idx = [0]*4; idx[axis] = 1; p.a[tuple(idx)] = c
        return p
    X = lin([1/np.sqrt(2*kappa), 1/np.sqrt(2*kappa), 0, 0])
    Z = lin([0, 0, 1/np.sqrt(2*nu), 1/np.sqrt(2*nu)])
    PR = lin([-1j*np.sqrt(kappa/2), 1j*np.sqrt(kappa/2), 0, 0])
    PZ = lin([0, 0, -1j*np.sqrt(nu/2), 1j*np.sqrt(nu/2)])
    # powers of X and Z
    Xp = [Poly4(N)]; Xp[0].a[0, 0, 0, 0] = 1
    Zp = [Poly4(N)]; Zp[0].a[0, 0, 0, 0] = 1
    for k in range(1, N+1):
        Xp.append(mul(Xp[-1], X)); Zp.append(mul(Zp[-1], Z))
    for (i, j), c in coefs.items():
        if i + j == 0 or i + j > N: continue
        H.a += c*mul(Xp[i], Zp[j]).a
    H.a += 0.5*(mul(PR, PR).a + mul(PZ, PZ).a)
    return H.truncate()

def birkhoff(H, kappa, nu, N):
    idx = np.indices(H.a.shape)
    a1, b1, a2, b2 = idx
    omega_dot = kappa*(a1 - b1) + nu*(a2 - b2)
    deg = idx.sum(axis=0)
    resonant = (a1 == b1) & (a2 == b2)
    chi_max = {}; nf_max = {}; chi_sum = {}
    Hc = H.copy()
    for n in range(3, N+1):
        Hn = Hc.degree_slice(n)
        chi = Poly4(N)
        mask = (deg == n) & (~resonant)
        with np.errstate(divide='ignore', invalid='ignore'):
            chi.a[mask] = 1j*Hn.a[mask]/omega_dot[mask]
        chi_max[n] = np.abs(chi.a).max(); chi_sum[n] = np.abs(chi.a).sum()
        # H <- exp(L_chi) H = sum_k {..{H,chi}..,chi}/k!
        term = Hc.copy(); Hnew = Hc.copy(); k = 1
        while True:
            term = bracket(term, chi); term.a /= k
            if np.abs(term.a).max() == 0 or k*(n-2) + 2 > N: 
                Hnew.a += term.a; break
            Hnew.a += term.a; k += 1
        Hc = Hnew.truncate()
        nf = Hc.degree_slice(n); nf_max[n] = np.abs(nf.a[resonant & (deg == n)]).max() if (resonant & (deg == n)).any() else 0
    return chi_max, chi_sum, nf_max

def run(label, Phi, N=14, R0=sp.Integer(1)):
    t0 = time.time()
    coefs, Lz2 = taylor_V(Phi, R0, N)
    kappa = np.sqrt(2*coefs[(2, 0)].real); nu = np.sqrt(2*coefs[(0, 2)].real)
    H = to_complex_vars(coefs, kappa, nu, N)
    chi_max, chi_sum, nf_max = birkhoff(H, kappa, nu, N)
    print(f"--- {label}: kappa={kappa:.5f} nu={nu:.5f} kappa/nu={kappa/nu:.5f}  (Taylor+BNF in {time.time()-t0:.0f}s)")
    prev = None
    for n in range(3, N+1):
        ratio = chi_max[n]/chi_max[n-1] if n > 3 and chi_max[n-1] > 0 else float('nan')
        print(f"   n={n:2d}  max|chi_n|={chi_max[n]:.3e}  ratio={ratio:8.3f}  ratio/n={ratio/n:7.3f}  max|NF_n|={nf_max[n]:.3e}")
    return chi_max

if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    r = sp.sqrt(R**2 + z**2)
    P2 = (3*z**2/r**2 - 1)/2
    cases = {
        'spherical Hernquist (integrable)': -1/(1+r),
        'Kuzmin-Kutuzov a=1 c=0.5 (integrable Staeckel)': -1/sp.sqrt(R**2 + z**2 + 1 + sp.Rational(1,4) + 2*sp.sqrt(sp.Rational(1,4) + sp.Rational(1,4)*R**2 + z**2)),
        'flattened Hernquist eps=0.3 (non-integrable)': -1/(1+r) + sp.Rational(3,10)*r**2/(1+r)**4*P2,
        'flattened core -1/(1+R^2+z^2/q^2), q=0.8': -1/(1 + R**2 + z**2/sp.Rational(16,25)),
    }
    which = sys.argv[2:] if len(sys.argv) > 2 else list(cases)
    for k in which:
        run(k, cases[k], N)
