import numpy as np, sympy as sp, sys
from ve3 import *
from ve3b import analyse2
A = -1/(1+r); B = r**2/(1+r)**4; E = -0.5
U, W, Y0 = axis_data(A, B, {2: 1}, 0)          # spherical Hernquist data
for label, dY in [('arbitrary dY = 0.3/(z(1+z)^3)', sp.Rational(3,10)/(z*(1+z)**3)), ('arbitrary dY = 0.5 z/(1+z)^6', sp.Rational(1,2)*z/(1+z)**6), ('arbitrary dY = 0.2/z^2', sp.Rational(1,5)/z**2)]:
    Y = sp.cancel(Y0 + dY)
    sing = singular_points(U, W, Y, E); Ms = monodromies(build_system(U, W, Y, E), sing, nseg=256, rtol=1e-13, atol=1e-15)
    analyse2(Ms, label); sys.stdout.flush()
