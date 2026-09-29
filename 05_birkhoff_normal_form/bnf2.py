import numpy as np, sys, time
from scipy.signal import fftconvolve
from bnf import Poly4, to_complex_vars, birkhoff

class S2:
    """truncated bivariate power series in (x,z), total degree <= N, dense (N+1)x(N+1)"""
    def __init__(self, N, a=None):
        self.N = N; self.a = np.zeros((N+1, N+1)) if a is None else a
    @staticmethod
    def const(N, c): s = S2(N); s.a[0, 0] = c; return s
    @staticmethod
    def var(N, which): s = S2(N); s.a[1, 0] = 1 if which == 'x' else 0; s.a[0, 1] = 1 if which == 'z' else 0; return s
    def trunc(self):
        i, j = np.indices(self.a.shape); self.a[i + j > self.N] = 0; return self
    def __add__(self, o): return S2(self.N, self.a + (o.a if isinstance(o, S2) else 0) + (0 if isinstance(o, S2) else np.eye(1, self.N+1, 0).T @ np.eye(1, self.N+1, 0)*0 + self._c(o)))
    def _c(self, c): m = np.zeros_like(self.a); m[0, 0] = c; return m
    def __radd__(self, o): return self.__add__(o)
    def __sub__(self, o): return self + (o*(-1) if isinstance(o, S2) else -o)
    def __rsub__(self, o): return (self*(-1)) + o
    def __mul__(self, o):
        if isinstance(o, S2):
            return S2(self.N, fftconvolve(self.a, o.a)[:self.N+1, :self.N+1]).trunc()
        return S2(self.N, self.a*o)
    __rmul__ = __mul__
    def inv(self):
        c0 = self.a[0, 0]; y = S2.const(self.N, 1/c0)
        for _ in range(int(np.ceil(np.log2(self.N+1))) + 1):
            y = y*(2 - self*y)
        return y
    def __truediv__(self, o): return self*o.inv() if isinstance(o, S2) else self*(1/o)
    def sqrt(self):
        c0 = self.a[0, 0]; y = S2.const(self.N, np.sqrt(c0))
        for _ in range(int(np.ceil(np.log2(self.N+1))) + 1):
            y = (y + self*y.inv())*0.5
        return y
    def __pow__(self, k):
        out = S2.const(self.N, 1.0)
        for _ in range(k): out = out*self
        return out

def build_V(case, N, R0=1.0, q=0.8, eps=0.3):
    x = S2.var(N, 'x'); z = S2.var(N, 'z')
    Rr = x + R0
    r2 = Rr*Rr + z*z
    r = r2.sqrt()
    if case == 'sph':
        Phi = lambda: (1 + r).inv()*(-1.0)
        dPhi_dR0 = 1/(1+R0)**2                        # Phi_R at (R0,0) for Phi=-1/(1+r)
    elif case == 'flat':
        P2 = (z*z*r2.inv()*3.0 - 1.0)*0.5
        Phi = lambda: (1 + r).inv()*(-1.0) + (r2*((1 + r)**4).inv()*P2)*eps
        # Phi_R at (R0,0): d/dR[-1/(1+R) + eps R^2/(1+R)^4 * (-1/2)]
        dPhi_dR0 = 1/(1+R0)**2 + eps*(-0.5)*(2*R0/(1+R0)**4 - 4*R0**2/(1+R0)**5)
    elif case == 'core':
        m2 = Rr*Rr + z*z*(1/q**2)
        Phi = lambda: (1 + m2).inv()*(-1.0)
        dPhi_dR0 = 2*R0/(1+R0**2)**2
    elif case == 'kk':
        a, c = 1.0, 0.5
        inner = (Rr*Rr*(c*c) + z*z*(a*a) + a*a*c*c).sqrt()
        Phi = lambda: (r2 + (a*a + c*c) + inner*2.0).sqrt().inv()*(-1.0)
        # numeric derivative at (R0,0)
        f = lambda RR: -1/np.sqrt(RR**2 + a*a + c*c + 2*np.sqrt(c*c*RR**2 + a*a*c*c))
        h = 1e-6; dPhi_dR0 = (f(R0+h) - f(R0-h))/(2*h)
    Ph = Phi()
    Lz2 = R0**3*Ph.a[1, 0]                      # equilibrium: -Lz^2/R0^3 + Phi_R = 0
    V = (Rr*Rr).inv()*(Lz2/2) + Ph
    coefs = {}
    for i in range(N+1):
        for j in range(N+1-i):
            if i + j > 0 and abs(V.a[i, j]) > 0: coefs[(i, j)] = complex(V.a[i, j])
    assert abs(coefs.get((1, 0), 0)) < 1e-5, ("not an equilibrium", coefs.get((1, 0)))
    coefs.pop((1, 0), None)
    return coefs

def run(case, label, N, lam=3.0):
    t0 = time.time()
    coefs = build_V(case, N)
    kappa = np.sqrt(2*coefs[(2, 0)].real); nu = np.sqrt(2*coefs[(0, 2)].real)
    H = to_complex_vars(coefs, kappa, nu, N)
    deg = np.indices(H.a.shape).sum(axis=0); H.a *= lam**(2.0 - deg)      # K(z) = lam^2 H(z/lam)
    chi_max, chi_sum, nf_max = birkhoff(H, kappa, nu, N)
    print(f"--- {label}: kappa={kappa:.5f} nu={nu:.5f} kappa/nu={kappa/nu:.5f} lam={lam} ({time.time()-t0:.0f}s)", flush=True)
    for n in range(4, N+1):
        ratio = chi_max[n]/chi_max[n-1]
        print(f"   n={n:2d}  max|chi_n|={chi_max[n]:.3e}  ratio={ratio:8.3f}  ratio/n={ratio/n:7.4f}   max|NF_n|={nf_max[n]:.3e}", flush=True)
    return chi_max, nf_max

if __name__ == '__main__':
    N = int(sys.argv[1]); lam = float(sys.argv[2]); cases = sys.argv[3:] or ['sph', 'flat', 'kk', 'core']
    labels = {'sph': 'spherical Hernquist (integrable)', 'flat': 'flattened Hernquist eps=0.3 (non-integrable)', 'kk': 'Kuzmin-Kutuzov (integrable)', 'core': 'flattened core q=0.8 (probe)'}
    for c in cases: run(c, labels[c], N, lam)
