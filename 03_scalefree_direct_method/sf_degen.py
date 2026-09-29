import sympy as sp
from sf_classify import classify, Fj
from sf_classify2 import canon, S, C
res = classify()
M = res['M']
det = M.det()
n, d = canon(det)
print("determinant (canonical numerator), factored:")
print(sp.factor(n))
