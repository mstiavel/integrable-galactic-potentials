import numpy as np, itertools, sympy as sp
from scipy.integrate import solve_ivp
zs = sp.symbols('z')

def family_F(c, eps):
    return (c*zs**2 - 1)/(1+zs), eps*zs**2/(1+zs)          # Hernquist cusp -> r^1 envelope
def family_G(c, eps):
    return (c*zs**3 - 1)/(1+zs**2), eps*zs**3/(1+zs**2)    # cored centre -> r^1 envelope

def nve(A, B, orbit):
    if orbit == 'axis': Us = sp.cancel(A + B); sgn = 1
    else:               Us = sp.cancel(A - B/2); sgn = -1
    dUs = sp.cancel(sp.diff(Us, zs))
    Ws = sp.cancel(dUs/zs - 3*sgn*B/zs**2)
    f = lambda e: sp.lambdify(zs, e, 'numpy')
    return f(Us), f(dUs), f(Ws), Us, Ws

def distinct_roots(expr_den):
    p = sp.Poly(expr_den, zs)
    if p.degree() <= 0: return []
    p = sp.Poly(sp.sqf_part(p.as_expr(), zs), zs)
    return [complex(r) for r in p.nroots(n=30)] if p.degree() > 0 else []

def singular_points(Us, Ws, E):
    num, den = sp.fraction(sp.cancel(E - Us))
    tps = distinct_roots(num)
    poles = distinct_roots(sp.fraction(Us)[1]) + distinct_roots(sp.fraction(Ws)[1])
    out = []
    for p in [0j] + tps + poles:
        if all(abs(p-q) > 1e-8 for q in out): out.append(p)
    return out

def integrate_path(U, dU, W, E, pts, rtol=1e-12, atol=1e-14):
    M = np.eye(2, dtype=complex)
    for z0, z1 in zip(pts[:-1], pts[1:]):
        dz = z1-z0
        def rhs(s, y):
            zz = z0+s*dz; Y = y.reshape(2, 2)
            den = 2*(E-U(zz)); a = dU(zz)/den; b = -W(zz)/den
            return (np.array([[0, 1], [b, a]]) @ Y * dz).reshape(-1)
        sol = solve_ivp(rhs, (0, 1), M.reshape(-1), method='DOP853', rtol=rtol, atol=atol)
        M = sol.y[:, -1].reshape(2, 2)
    return M

def sl2(M): return M/np.sqrt(np.linalg.det(M))
def comm(P_, Q): return P_ @ Q @ np.linalg.inv(P_) @ np.linalg.inv(Q)

def analyse(A, B, E, orbit, base=0.6+0.9j, frac=0.4, nseg=96):
    U, dU, W, Us, Ws = nve(A, B, orbit)
    sing = singular_points(Us, Ws, E)
    order = np.argsort([np.angle(s-base) for s in sing])
    Ms = []
    for k in order:
        s = sing[k]
        others = [t for t in sing if abs(t-s) > 1e-9]
        rad = frac*min([abs(s-t) for t in others] + [abs(s-base)])
        d = (base-s)/abs(base-s); start = s+rad*d; th0 = np.angle(d)
        circ = [s+rad*np.exp(1j*(th0+2*np.pi*j/nseg)) for j in range(nseg+1)]
        Ms.append(integrate_path(U, dU, W, E, [base, start]+circ[1:]+[start, base]))
    Pinf = np.eye(2, dtype=complex)
    for M in Ms: Pinf = M @ Pinf
    S = [sl2(M) for M in Ms]
    comms = [comm(S[i], S[j]) for i, j in itertools.combinations(range(len(S)), 2)]
    relC = 0.0 if len(comms) < 2 else max(np.linalg.norm(C1@C2-C2@C1)/(np.linalg.norm(C1)*np.linalg.norm(C2)) for C1, C2 in itertools.combinations(comms, 2))
    hyper = max(abs(np.trace(S[i]@S[j])) for i, j in itertools.permutations(range(len(S)), 2))
    return dict(sing=sing, Pinf_eig=np.linalg.eigvals(Pinf), derived=relC, hyper=hyper, maxnorm=max(np.linalg.norm(M) for M in Ms))

def predicted_inf_eigs(mu, a, b, orbit):
    u, be = (a+b, b) if orbit == 'axis' else (a-b/2, -b)
    D = np.sqrt(complex((mu+2)**2 - 24*be/u))
    rho = [((2-mu)+D)/4, ((2-mu)-D)/4]
    return np.exp(-2j*np.pi*np.array(rho))   # loop around infinity traversed as product of finite loops

if __name__ == '__main__':
    c = 0.5
    for fam, name in [(family_F, 'F: Hernquist cusp -> r^1 envelope'), (family_G, 'G: cored centre -> r^1 envelope')]:
        print("\n######", name)
        for eps in [0.0, 0.25, -0.2, 0.15, 0.35, -0.1, -0.3]:
            for orbit in ['axis', 'eq']:
                A, B = fam(c, eps)
                for E in [0.5, 2.0]:
                    r = analyse(A, B, E, orbit)
                    pe = predicted_inf_eigs(1.0, c, eps, orbit)
                    got = np.sort_complex(np.round(r['Pinf_eig'], 4)); exp_ = np.sort_complex(np.round(pe, 4))
                    print(f" eps={eps:+.2f} {orbit:4s} E={E:.1f} | loop@inf eig {got} predicted {exp_} | derived-subgroup noncomm {r['derived']:.2e} | max|tr| {r['hyper']:.2f} | max|M| {r['maxnorm']:.1f}")
