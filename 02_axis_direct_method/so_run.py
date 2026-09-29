import sys, sympy as sp
from second_order import run, b1_family_dim
N, M = int(sys.argv[1]), int(sys.argv[2])
print(f"=== alpha=1, beta=0, n4, N={N}, M={M}", flush=True)
conds, tpar, lowb, sub1, b1 = run(N, M, 1, 0)
for d, vec in b1_family_dim(conds, tpar, lowb, b1, M):
    print(f"  surviving B1 family: dimension {d}; b2,b4,b6,b8,b10,b12 =", [vec[i] for i in [1,3,5,7,9,11] if i < len(vec)], flush=True)
