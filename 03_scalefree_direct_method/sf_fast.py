import sympy as sp, sys, time, itertools, pickle, random
from sf_general import build, Fj, J, S, C, redSC, th

def to_SC(e):
    """expression in sin(th), cos(th) (and rational) -> (numerator poly, monomial denominator) in S, C"""
    e = sp.together(sp.expand(e.subs({sp.sin(th): S, sp.cos(th): C})))
    n, d = sp.fraction(e)
    return redSC(n), sp.factor(d)

def Dsc(e):
    out = sp.diff(e, S)*C - sp.diff(e, C)*S
    return out

class Elim:
    def __init__(self, n, d, muval, verbose=True, kunit=None):
        t0 = time.time()
        nparams = len([1 for i in range(n-d, -1, -1) if i % 2 == 0 and (n-d-i) % 2 == 0])
        self.nparams = nparams
        params = sp.symbols(f'a0:{nparams}')
        if kunit is not None: params = tuple(sp.Integer(1) if i == kunit else sp.Integer(0) for i in range(nparams))
        self.params = params
        eqsJ, Uj, unknowns = build(n, d, params, muval)
        self.Uj = Uj
        # convert equations to S,C polynomials (they are polynomial in sin,cos after multiplying by r-powers)
        eqs = {}
        for key, e in eqsJ.items():
            num, den = to_SC(e); eqs[key] = sp.expand(num)   # den is 1 or monomial; multiply through
        self.eqs = eqs
        if verbose: print(f"  equations in (S,C) in {time.time()-t0:.0f}s", flush=True)
        jetsyms = [Uj[key][k] for key in Uj for k in range(J)]
        def D(e):
            out = Dsc(e)
            for key in Uj:
                for k in range(J-1): out += sp.diff(e, Uj[key][k])*Uj[key][k+1]
            for k in range(J-1): out += sp.diff(e, Fj[k])*Fj[k+1]
            return sp.expand(out)
        self.D = D
        # level-by-level: express first derivatives of unknowns; keep (num, den) with monomial den
        subsP = {}
        rels = []
        levels = sorted({key[0] for key in eqs}, reverse=True)
        def apply(e):
            # substitute u' -> num/den  and clear: return expanded numerator over common monomial denominator
            e = sp.expand(e.subs(subsP))
            e = sp.together(e); num, den = sp.fraction(e)
            return redSC(sp.expand(num)), den
        for m in levels:
            if m == 1:
                D0, D1 = Uj[(0, 0)][0], Uj[(0, 0)][1]
                E1, _ = apply(eqs[(1, 1, 0)]); E0, _ = apply(eqs[(1, 0, 1)])
                # E1 = C*D' + wf*S*D - T1 ; E0 = -S*D' + wf*C*D - T0  (up to overall monomial factors); extract wf
                cD1 = sp.Poly(E1, D1).coeff_monomial(D1); cD0 = sp.Poly(E1, D0).coeff_monomial(D0)
                wf = sp.cancel(cD0/(cD1*S/C))
                T1 = -sp.expand(E1.subs({D0: 0, D1: 0}))/cD1*C; T0e = sp.expand(E0.subs({D0: 0, D1: 0}))
                cD1b = sp.Poly(E0, D1).coeff_monomial(D1)
                T0 = -T0e/cD1b*(-S)
                Dsol = sp.cancel((S*T1 + C*T0)/wf); Dp = sp.cancel(C*T1 - S*T0)
                rel = sp.together(D(Dsol) - Dp); rel = sp.expand(rel.subs(subsP)); rel = sp.together(rel)
                rels.append(redSC(sp.numer(rel)))
                continue
            for a in range(m, 0, -1):
                e, _ = apply(eqs[(m, a, m-a)])
                var = Uj[(m-1, a-1)][1]
                co = sp.Poly(e, var).coeff_monomial(var); rest = sp.expand(e - co*var)
                subsP[var] = sp.cancel(-rest/co)
                for k2 in list(subsP): subsP[k2] = sp.cancel(subsP[k2].subs(var, subsP[var]))
            e, _ = apply(eqs[(m, 0, m)])
            rels.append(e)
        self.subsP = subsP; self.rels = rels
        if verbose: print(f"  {len(rels)} relations in {time.time()-t0:.0f}s; sizes {[len(str(r)) for r in rels]}", flush=True)
        self.u0 = [Uj[key][0] for key in Uj if key[0] > 0]
        nu = len(self.u0)
        # derivative pool
        def red(e):
            e = sp.together(sp.expand(e.subs(subsP))); return redSC(sp.numer(e))
        pool = []; tags = []
        counts = [nu - 2*(len(rels) - 1) + 1] + [4]*(len(rels) - 1)
        for i, (R, cnt) in enumerate(zip(rels, counts)):
            cur = red(R)
            for k in range(cnt):
                pool.append(cur); tags.append((i, k)); cur = red(D(cur))
        M, rhs = sp.linear_eq_to_matrix(pool, self.u0)
        self.M, self.rhs, self.tags = M, rhs, tags
        if verbose: print(f"  jet matrix {M.shape} in {time.time()-t0:.0f}s", flush=True)
        # choose rows and columns numerically (matrix may be rank deficient: free homogeneous solutions)
        rng = random.Random(3)
        num = {s_: sp.Rational(rng.randint(2, 40), rng.randint(2, 40)) for s_ in list(Fj) + [S, C]}
        Mn = M.subs(num)
        rk = Mn.rank()
        if verbose: print(f"  numeric rank {rk} of {nu} unknowns", flush=True)
        # independent columns
        cols = []
        for j in range(nu):
            if Mn.extract(list(range(M.shape[0])), cols + [j]).rank() > len(cols): cols.append(j)
        rows = []
        for i in range(M.shape[0]):
            if Mn.extract(rows + [i], cols).rank() > len(rows): rows.append(i)
        self.best = tuple(rows); self.cols = cols
        free = [self.u0[j] for j in range(nu) if j not in cols]
        if verbose: print("  solve rows:", [tags[i] for i in rows], " free unknowns:", free, flush=True)
        from sympy.polys.matrices import DomainMatrix
        sub = M.extract(rows, cols)
        rsub = rhs.extract(rows, [0]) - M.extract(rows, [j for j in range(nu) if j not in cols])*sp.Matrix(free) if free else rhs.extract(rows, [0])
        allsyms = set()
        for e in list(sub) + list(rsub) + list(M.extract([i for i in range(M.shape[0]) if i not in rows], list(range(nu)))) + list(rhs):
            allsyms |= e.free_symbols
        gens = sorted(allsyms, key=str)
        if verbose: print("  field generators:", gens, flush=True)
        K_ = sp.QQ.frac_field(*gens)
        A = DomainMatrix.from_Matrix(sub).convert_to(K_); B = DomainMatrix.from_Matrix(rsub).convert_to(K_)
        X = A.lu_solve(B)
        self.rest = [i for i in range(M.shape[0]) if i not in rows]
        if verbose: print(f"  solved in {time.time()-t0:.0f}s", flush=True)
        # residual conditions computed inside the field domain
        Mrest = DomainMatrix.from_Matrix(M.extract(self.rest, cols)).convert_to(K_)
        rrest = rhs.extract(self.rest, [0]) - (M.extract(self.rest, [j for j in range(nu) if j not in cols])*sp.Matrix(free) if free else sp.zeros(len(self.rest), 1))
        Rr = DomainMatrix.from_Matrix(rrest).convert_to(K_)
        Res = (Mrest*X - Rr).to_Matrix()
        conds = []
        for v in Res:
            num = sp.numer(v)
            num = sp.factor(redSC(sp.expand(num)))
            if any(num.has(fu) for fu in free): print("   WARNING: condition depends on free unknown", flush=True)
            conds.append(num)
        self.conds = conds
        if verbose:
            print(f"  {len(conds)} conditions in {time.time()-t0:.0f}s", flush=True)
            for c in conds: print("   COND:", str(c)[:1200], flush=True)

if __name__ == '__main__':
    n, d, muval = int(sys.argv[1]), int(sys.argv[2]), sp.Rational(sys.argv[3])
    kunit = int(sys.argv[4]) if len(sys.argv) > 4 else None
    E = Elim(n, d, muval, kunit=kunit)
    pickle.dump((E.conds, E.params), open(f'fast_deg{n}_d{d}_mu{str(muval).replace("/","o")}' + (f'_k{kunit}' if kunit is not None else '') + '.pkl', 'wb'))
