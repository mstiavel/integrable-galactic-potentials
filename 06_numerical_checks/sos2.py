import numpy as np, sys
from sos import run
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ts = [float(a) for a in sys.argv[1:]]
fig, axs = plt.subplots(1, len(ts), figsize=(6*len(ts), 5.5))
axs = np.atleast_1d(axs)
for ax, t in zip(axs, ts):
    lyap, pts = run(t, norb=28, T=600.0)
    for xs, pxs in pts: ax.plot(xs, pxs, '.', ms=1.2)
    ax.set_title(f"t={t}: median={np.median(lyap):.4f}, max={lyap.max():.4f}, n(>0.02)={int((lyap>0.02).sum())}")
    ax.set_xlabel('x (z=0, pz>0)'); ax.set_ylabel('px')
    print(f"t={t}: lyap sorted: {np.round(np.sort(lyap),4)}", flush=True)
plt.tight_layout(); plt.savefig(f'/mnt/user-data/outputs/sos_{"_".join(str(t) for t in ts)}.png', dpi=100)
