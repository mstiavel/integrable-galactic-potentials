import sympy as sp, sys
from n6_jet import eliminate_jets, indicial
z = sp.symbols('z'); a, b, d = sp.symbols('alpha beta delta')
cases = [("n=6, A=z, degenerate seed alpha=0 (beta,delta)", 6, [0, b, d], z, 5),
         ("n=6, Jaffe local A=ln z, seed alpha (L^2-type)", 6, [a], sp.log(z), 5),
         ("n=4, Jaffe local A=ln z, seed alpha (check vs earlier (s+1)(s+2)^2)", 4, [a], sp.log(z), 4)]
for label, n, seed, A, D in cases:
    res = eliminate_jets(n, seed, A, D)
    if res is None: print(label, ": no elimination"); continue
    bexprs, bj = res
    order, ind = indicial(bexprs[0], bj)
    print(f"{label}: ODE order {order}, indicial: {ind}")
    sys.stdout.flush()
