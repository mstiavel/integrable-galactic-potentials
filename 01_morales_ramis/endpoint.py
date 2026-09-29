"""
Endpoint condition.  At a scale-free endpoint  A ~ a r^mu,  B ~ b r^mu  (potential-dominated,
i.e. r->inf for mu>0), the NVE exponents differ by  delta = Delta/2,
     Delta^2 = (mu+2)^2 - 24 b_eff/u_eff,
with (u_eff, b_eff) = (a+b, b) on the axis and (a-b/2, -b) in the plane.
Morales-Ramis admissible values of delta for homogeneous potentials of degree mu (general-mu
families of the MR table; specific rational mu admit finitely many extra exceptional values):
     delta = mu (p + 1/2)                 [family  (k, 1/2((k-1)/k + p(p+1)k))]
     delta = | mu (p - 1/2) + 1 |         [family  (k, p + p(p-1)k/2)]
Spherical case b=0: Delta = mu+2, delta = mu/2+1 = family 2 with p=1  (consistent).
Also: lambda_MR = 1 - 3 b_eff/(mu u_eff), and Delta^2 = (mu-2)^2 + 8 mu lambda  (Yoshida form).
"""
import numpy as np
from scipy.optimize import brentq

def admissible(mu, pmax=6):
    vals = []
    for p in range(-pmax, pmax+1):
        vals.append(('T', p, mu*(p+0.5)))          # family 18-type
        vals.append(('L', p, abs(mu*(p-0.5)+1)))   # family 1-type
    return [v for v in vals if v[2] >= 0]

def t_from_delta_axis(mu, d):
    # 24 t/(1+t) = (mu+2)^2 - 4 d^2  =: s   ->  t = s/(24-s)
    s = (mu+2)**2 - 4*d**2
    if abs(24-s) < 1e-12: return None
    return s/(24-s)

def delta_eq(mu, t):
    D2 = (mu+2)**2 + 24*t/(1-t/2)
    return np.sqrt(D2) if D2 >= 0 else None

def density_ok(mu, t):
    # Laplacian of r^mu (1 + t P2):  r^(mu-2) [ mu(mu+1) + t (mu(mu+1)-6) P2 ] >= 0 for P2 in [-1/2, 1]
    m = mu*(mu+1)
    return min(m + t*(m-6), m - 0.5*t*(m-6)) >= -1e-12

sols = []
mus = np.linspace(0.02, 2.0, 1981)
for fa, pa, _ in admissible(1.0):
    for fe, pe, _ in admissible(1.0):
        def g(mu):
            da = dict((f, p) for f, p, _ in [(fa, pa, 0)])
            dax = mu*(pa+0.5) if fa == 'T' else abs(mu*(pa-0.5)+1)
            t = t_from_delta_axis(mu, dax)
            if t is None or t <= -1 or t >= 1.99: return np.nan
            de = delta_eq(mu, t)
            if de is None: return np.nan
            dadm = mu*(pe+0.5) if fe == 'T' else abs(mu*(pe-0.5)+1)
            return de/2 - dadm
        vals = np.array([g(m) for m in mus])
        for i in range(len(mus)-1):
            if np.isnan(vals[i]) or np.isnan(vals[i+1]): continue
            if vals[i] == 0 or vals[i]*vals[i+1] < 0:
                try:
                    mu0 = brentq(g, mus[i], mus[i+1])
                except Exception:
                    continue
                dax = mu0*(pa+0.5) if fa == 'T' else abs(mu0*(pa-0.5)+1)
                t = t_from_delta_axis(mu0, dax)
                sols.append((round(mu0, 6), round(t, 6), fa, pa, fe, pe, density_ok(mu0, t)))
# also count exact identities g==0 on whole intervals (spherical t=0 lines) separately
uniq = {}
for s in sols:
    key = (s[0], s[1])
    uniq.setdefault(key, s)
print("mu (=2-gamma_out)   t=b/a     axis-family(p)  eq-family(q)  positive density?")
for k in sorted(uniq):
    s = uniq[k]
    if abs(s[1]) < 1e-9: continue   # spherical
    if abs(s[1]) > 1.0: continue
    print(f"{s[0]:8.4f}   {s[1]:9.5f}     {s[2]}({s[3]:+d})          {s[4]}({s[5]:+d})        {s[6]}")

# spot check mu=1 analytic: t=1/2 and t=-2/5
for t in [0.5, -0.4]:
    dax = np.sqrt(9 - 24*t/(1+t))/2; de = delta_eq(1.0, t)/2
    print(f"mu=1, t={t}: delta_axis={dax:.4f} delta_eq={de:.4f}  density_ok={density_ok(1.0,t)}")
# thread-2-like slope
for mu in [5/6, 1.0, 0.5, 1.5]:
    print(f"mu={mu:.4f}: admissible deltas (>=0, small p):", sorted(set(round(v,4) for _,_,v in admissible(mu, 3))))
