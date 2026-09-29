import numpy as np, sys
from scipy.integrate import solve_ivp
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

def make(t):
    def acc(x, z):
        r = np.hypot(x, z) + 1e-300
        c = z/r; P2 = (3*c*c-1)/2
        # Phi = r (1 + t P2);  P2 = (3 z^2/r^2 - 1)/2 -> Phi = r(1 - t/2) + 1.5 t z^2 / r
        dPhidx = x/r*(1 - t/2) - 1.5*t*z*z*x/r**3
        dPhidz = z/r*(1 - t/2) + 1.5*t*(2*z/r - z**3/r**3)
        return -dPhidx, -dPhidz
    def pot(x, z):
        r = np.hypot(x, z); return r*(1 - t/2) + (1.5*t*z*z/r if r > 0 else 0.0)
    return acc, pot

def rhs_factory(acc):
    def rhs(s, y):
        x, z, px, pz, dx, dz, dpx, dpz = y
        ax, az = acc(x, z)
        # tangent: numerical Jacobian of acc (central differences, cheap in 2D)
        h = 1e-6
        axp, azp = acc(x+h, z); axm, azm = acc(x-h, z)
        axq, azq = acc(x, z+h); axr, azr = acc(x, z-h)
        J = np.array([[(axp-axm)/(2*h), (axq-axr)/(2*h)], [(azp-azm)/(2*h), (azq-azr)/(2*h)]])
        d2 = J @ np.array([dx, dz])
        return [px, pz, ax, az, dpx, dpz, d2[0], d2[1]]
    return rhs

def run(t, E=1.0, norb=36, T=600.0, seed=1):
    acc, pot = make(t); rhs = rhs_factory(acc)
    rng = np.random.default_rng(seed)
    lyap = []; pts = []
    for k in range(norb):
        # launch from z=0 plane, pz>0, x in (0.1, xmax), px chosen from energy with fraction of radial energy
        while True:
            x0 = rng.uniform(0.15, 1.6); f = rng.uniform(0, 1)
            V = pot(x0, 0.0)
            if E - V > 0: break
        K = E - V; px0 = np.sqrt(2*K*f)*rng.choice([-1, 1]); pz0 = np.sqrt(2*K*(1-f))
        d = rng.normal(size=4); d /= np.linalg.norm(d)
        y = [x0, 0.0, px0, pz0, *d]
        tt = 0.0; ssum = 0.0; crossings = []
        def event(s, y): return y[1]
        event.direction = 1
        sol = solve_ivp(rhs, (0, T), y, method='DOP853', rtol=1e-10, atol=1e-12, events=event, dense_output=False, max_step=0.05)
        # Lyapunov from final tangent growth (single renormalisation is crude but fine for classification)
        dfin = np.array(sol.y[4:, -1]); growth = np.log(np.linalg.norm(dfin))/T
        lyap.append(growth)
        ev = sol.y_events[0]
        pts.append((ev[:, 0], ev[:, 2]))
    return np.array(lyap), pts

if __name__ == '__main__':
    ts = [0.0, -0.4, 0.5, 0.25]
    fig, axs = plt.subplots(1, 4, figsize=(20, 5))
    for ax, t in zip(axs, ts):
        lyap, pts = run(t)
        for xs, pxs in pts: ax.plot(xs, pxs, ',', ms=1)
        ax.set_title(f"t={t}: median lyap={np.median(lyap):.4f}, max={lyap.max():.4f}")
        ax.set_xlabel('x (z=0, pz>0)'); ax.set_ylabel('px')
        print(f"t={t}: finite-time Lyapunov (per unit time), sorted: {np.round(np.sort(lyap),4)}", flush=True)
    plt.tight_layout(); plt.savefig('/mnt/user-data/outputs/sos_scalefree.png', dpi=110)
