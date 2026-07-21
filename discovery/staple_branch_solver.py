"""High-precision staple branch curves over u = tan(beta/2), for the interval
theorem. Newton at 30 digits on the active-set systems:

  branch1 (u <= kink ~ base 39.05deg):  D_s = D_p3 = 0, stationarity det = 0
  vertex  (near kink):                  D_s = D_p3 = D_p4 = 0
  branch2 (u >= kink):                  D_s = D_p4 = 0, stationarity det = 0

where D_s is the self-symmetric near-active edge and D_p3 / D_p4 the two
mirror pairs. Produces per-segment polynomial fits with rational
coefficients, scaled by (1+eps) (scaling H about the origin enlarges it,
so every escape margin can only gain), and float sanity margins.
"""
import itertools, json
import numpy as np
import mpmath as mp
from scipy.spatial import ConvexHull

mp.mp.dps = 35

# ------------------------------------------------------------------ numerics
def tri_mp(u):
    u = mp.mpf(u)
    den = 1 + u*u
    sb = 2*u/den
    cb = (1 - u*u)/den
    return sb, cb

def rotmat(n):
    return [[n[0], n[1]], [-n[1], n[0]]]

def W_pts_mp(u, a, b, c):
    sb, cb = tri_mp(u)
    H = [(-b, mp.mpf(0)), (b, mp.mpf(0)), (a, c), (-a, c)]
    nws = [((mp.mpf(0), mp.mpf(-1)), 2*cb), ((sb, cb), mp.mpf(1)), ((-sb, cb), mp.mpf(1))]
    copies = []
    for n, w in nws:
        R = rotmat(n)
        copies.append([(w*(R[0][0]*v[0] + R[0][1]*v[1]), w*(R[1][0]*v[0] + R[1][1]*v[1]))
                       for v in H])
    pts = [ (p[0]+q[0]+r[0], p[1]+q[1]+r[1])
            for p, q, r in itertools.product(*copies) ]
    rho = 2*sb*cb
    return pts, rho

def cycle_float(u, a, b, c):
    pts, _ = W_pts_mp(u, a, b, c)
    P = np.array([[float(x), float(y)] for x, y in pts])
    return [int(i) for i in ConvexHull(P).vertices]

def margins_mp(u, a, b, c, cyc):
    pts, rho = W_pts_mp(u, a, b, c)
    Q = [pts[i] for i in cyc]
    n = len(Q)
    out = []
    for k in range(n):
        p, q = Q[k], Q[(k+1) % n]
        e = (q[0]-p[0], q[1]-p[1])
        cr = p[0]*q[1] - p[1]*q[0]
        out.append(cr*cr - rho*rho*(e[0]**2 + e[1]**2))
    return out

def length_mp(a, b, c):
    return 2*b + 2*mp.sqrt((a-b)**2 + c**2)

def grad(f, x, h=mp.mpf(10)**-12):
    g = []
    for i in range(3):
        xp = list(x); xm = list(x)
        xp[i] += h; xm[i] -= h
        g.append((f(xp) - f(xm))/(2*h))
    return g

def det3(A):
    return (A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
          - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
          + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]))

def solve_branch(u, x0, cyc, idx_s, idx_p, mode):
    """mode 'stat': tangency at idx_s, idx_p + stationarity; mode 'vertex':
    tangency at idx_s, idx_p[0], idx_p[1] (three distinct constraints)."""
    u = mp.mpf(u)
    def Dk(x, k):
        return margins_mp(u, x[0], x[1], x[2], cyc)[k]
    def F(x):
        if mode == 'stat':
            gL = grad(lambda y: length_mp(*y), x)
            g1 = grad(lambda y: Dk(y, idx_s), x)
            g2 = grad(lambda y: Dk(y, idx_p[0]), x)
            return [Dk(x, idx_s), Dk(x, idx_p[0]), det3([gL, g1, g2])]
        else:
            return [Dk(x, idx_s), Dk(x, idx_p[0]), Dk(x, idx_p[1])]
    x = [mp.mpf(v) for v in x0]
    for _ in range(60):
        Fx = F(x)
        J = mp.matrix(3, 3)
        h = mp.mpf(10)**-12
        for j in range(3):
            xp = list(x); xm = list(x)
            xp[j] += h; xm[j] -= h
            Fp, Fm = F(xp), F(xm)
            for i in range(3):
                J[i, j] = (Fp[i] - Fm[i])/(2*h)
        dx = mp.lu_solve(J, mp.matrix([-v for v in Fx]))
        x = [x[i] + dx[i] for i in range(3)]
        if max(abs(v) for v in dx) < mp.mpf(10)**-28:
            break
    return x

def classify_actives(u, x, cyc):
    mg = [float(v) for v in margins_mp(mp.mpf(u), *[mp.mpf(v) for v in x], cyc)]
    order = np.argsort(mg)
    # single = self-mirror edge (unique margin), pairs share values
    singles, pairs = [], []
    used = set()
    for i in range(12):
        if i in used: continue
        twins = [j for j in range(12) if j != i and abs(mg[j]-mg[i]) < 1e-12]
        if twins:
            pairs.append((i, twins[0], mg[i])); used.update({i, twins[0]})
        else:
            singles.append((i, mg[i])); used.add(i)
    singles.sort(key=lambda t: t[1]); pairs.sort(key=lambda t: t[2])
    return singles, pairs, mg

if __name__ == '__main__':
    rows = json.load(open('results/staple_fit_data.json'))
    by_u = {round(r['u'], 6): r for r in rows}

    U1, US1, US2, U2 = 0.278, 0.352, 0.356, 0.3885

    def seed_for(u):
        rr = min(rows, key=lambda r: abs(r['u'] - u))
        return [rr['a'], rr['b'], rr['c']]

    out = {}
    for name, ulo, uhi, mode, n in [
            ('branch1', U1 - 0.002, US1 + 0.0005, 'stat', 40),
            ('vertex',  US1 - 0.0005, US2 + 0.0005, 'vertex', 9),
            ('branch2', US2 - 0.0005, U2 + 0.002, 'stat', 40)]:
        grid = np.linspace(ulo, uhi, n)
        recs = []
        x_prev = None
        for u in grid:
            x0 = x_prev if x_prev is not None else seed_for(u)
            cyc = cycle_float(u, *[float(v) for v in x0])
            singles, prs, mg = classify_actives(u, x0, cyc)
            idx_s = singles[0][0]
            if mode == 'stat':
                idx_p = (prs[0][0], prs[0][1])
            else:
                idx_p = (prs[0][0], prs[1][0])
            x = solve_branch(u, x0, cyc, idx_s, idx_p, mode)
            recs.append(dict(u=float(u),
                             a=mp.nstr(x[0], 30), b=mp.nstr(x[1], 30),
                             c=mp.nstr(x[2], 30),
                             L=mp.nstr(length_mp(*x), 30),
                             idx=(idx_s, idx_p)))
            x_prev = [float(v) for v in x]
        out[name] = recs
        print(f"{name}: {len(recs)} pts, u in [{grid[0]:.4f},{grid[-1]:.4f}], "
              f"L from {recs[0]['L'][:12]} to {recs[-1]['L'][:12]}")
    json.dump(out, open('results/branch_curves.json', 'w'), indent=1)
