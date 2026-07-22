"""Diagonal-type (zigzag/'zee') vs boundary-type (staple/line-arc) escape
branches for isosceles triangles (legs 1, base angle beta): sweeps, the
regime crossing near 42.3 deg, and its transversality (slope gap).

Diagonal path: A -> tangent -> arc(C1,r) -> inner tangent -> arc(C2,r)
-> tangent -> D;  r = 0 is the polygonal 3-segment zigzag (Besicovitch zee).
All optimizations use cutting planes: SLSQP on a working direction grid,
then dense-grid violations are appended until the dense margin is clean.
Small residual deficits eps are repaired exactly by scaling the hull by
rho/(rho-eps) (W is positively 1-homogeneous), correcting L accordingly.
"""
import json
import warnings
import numpy as np
from scipy.optimize import minimize
import zalgalloid_family as ZF

np.seterr(all='ignore')
warnings.filterwarnings('ignore')

ZETA = 2.2782916414368529
THD = np.linspace(0, 2*np.pi, 200001)[:-1]

def u_of(th):
    return np.array([np.cos(th), np.sin(th)])

def cross2(a, b):
    return a[0]*b[1] - a[1]*b[0]

def rho_of(beta):
    return 2*np.sin(beta)*np.cos(beta)

def rot_mats(beta):
    sb, cb = np.sin(beta), np.cos(beta)
    return (np.array([[0, -1], [1, 0]]),
            np.array([[sb, cb], [-cb, sb]]),
            np.array([[-sb, cb], [-cb, -sb]]), sb, cb)

# ---------------- diagonal (zigzag) branch --------------------------------

def zig_geom(p):
    """p = [Ax,Ay,C1x,C1y,C2x,C2y,Dx,Dy,r] -> features dict or None."""
    A = p[0:2]; C1 = p[2:4]; C2 = p[4:6]; D = p[6:8]; r = p[8]
    if r < 1e-10:
        pts = np.array([A, C1, C2, D])
        L = (np.linalg.norm(C1-A) + np.linalg.norm(C2-C1) + np.linalg.norm(D-C2))
        return dict(r=0.0, pts=pts, L=L, arcs=[])
    dA = np.linalg.norm(A - C1); dD = np.linalg.norm(D - C2)
    dm = np.linalg.norm(C2 - C1)
    if dA <= r*1.0001 or dD <= r*1.0001 or dm <= 2*r*1.0001:
        return None
    s1 = np.sign(cross2(C1 - A, C2 - C1))
    if s1 == 0:
        return None
    s2 = -s1
    phiA = np.arctan2(*(A - C1)[::-1]); aA = np.arccos(r / dA)
    phiD = np.arctan2(*(D - C2)[::-1]); aD = np.arccos(r / dD)
    psi = np.arctan2(*(C2 - C1)[::-1]); am = np.arccos(2*r / dm)
    def pick(cands, C, ref, s, entry):
        best = None
        for th in cands:
            T = C + r*u_of(th)
            t = s*np.array([-np.sin(th), np.cos(th)])
            v = (T - ref) if entry else (ref - T)
            nv = np.linalg.norm(v)
            if nv < 1e-14:
                continue
            err = abs(cross2(v/nv, t))
            if v @ t > 0 and (best is None or err < best[0]):
                best = (err, th, T)
        return best
    mid = None
    for sgn in (1, -1):
        th1 = psi + sgn*am; th2 = th1 + np.pi
        T1 = C1 + r*u_of(th1); T2 = C2 + r*u_of(th2)
        v = T2 - T1; nv = np.linalg.norm(v)
        if nv < 1e-12:
            continue
        t1 = s1*np.array([-np.sin(th1), np.cos(th1)])
        t2 = s2*np.array([-np.sin(th2), np.cos(th2)])
        if v/nv @ t1 > 0.999999 and v/nv @ t2 > 0.999999:
            mid = (th1, th2, T1, T2, nv)
            break
    ein = pick([phiA + aA, phiA - aA], C1, A, s1, True)
    eout = pick([phiD + aD, phiD - aD], C2, D, s2, False)
    if mid is None or ein is None or eout is None or ein[0] > 1e-8 or eout[0] > 1e-8:
        return None
    th_in1, T_in1 = ein[1], ein[2]
    th_out1, th_in2, T_out1, T_in2, tm = mid
    th_out2, T_out2 = eout[1], eout[2]
    sp1 = np.mod(s1*(th_out1 - th_in1), 2*np.pi)
    sp2 = np.mod(s2*(th_out2 - th_in2), 2*np.pi)
    if sp1 > np.pi or sp2 > np.pi:
        return None
    L = (np.sqrt(dA**2 - r**2) + r*sp1 + tm + r*sp2 + np.sqrt(dD**2 - r**2))
    pts = np.array([A, T_in1, T_out1, T_in2, T_out2, D])
    arcs = [(C1, th_in1, sp1, s1), (C2, th_in2, sp2, s2)]
    return dict(r=r, pts=pts, L=L, arcs=arcs)

def support_zig(g, U):
    h = (U @ g['pts'].T).max(axis=1)
    r = g['r']
    if r > 0:
        ang = np.arctan2(U[:, 1], U[:, 0])
        for C, th0, span, s in g['arcs']:
            rel = np.mod(s*(ang - th0), 2*np.pi)
            m = rel <= span
            h = np.where(m, np.maximum(h, U @ C + r), h)
    return h

def margins_zig(beta, g, thetas):
    R0, R1, R2, sb, cb = rot_mats(beta)
    U = np.stack([np.cos(thetas), np.sin(thetas)], 1)
    return (2*cb*support_zig(g, U @ R0) + support_zig(g, U @ R1)
            + support_zig(g, U @ R2)) - 2*sb*cb

def zee_classic(beta):
    sb, cb = np.sin(beta), np.cos(beta)
    L = 6*sb*cb/np.sqrt(1 + 8*sb**2)
    delta = np.arctan(1/(3*np.tan(beta)))
    pts = [np.zeros(2)]
    for sgn in (1, -1, 1):
        pts.append(pts[-1] + (L/3)*np.array([np.cos(sgn*delta), np.sin(sgn*delta)]))
    return L, np.array(pts)

def cutting_opt(beta, p0, cost, geom, marg, bounds, rounds=6, tol=1e-7):
    """generic cutting-plane SLSQP; marg(p, thetas) -> margins array."""
    TH = np.linspace(0, 2*np.pi, 4001)[:-1]
    p = np.asarray(p0, float)
    for _ in range(rounds):
        res = minimize(cost, p, constraints=[dict(type='ineq',
                       fun=lambda q: marg(q, TH))],
                       bounds=bounds, method='SLSQP',
                       options=dict(ftol=1e-14, maxiter=600))
        p = res.x
        m = marg(p, THD)
        mm = m.min()
        if mm > -tol:
            break
        bad = THD[m < max(mm*0.3, -tol)]
        if len(bad) > 400:
            bad = bad[::len(bad)//400 + 1]
        TH = np.sort(np.r_[TH, bad])
    mm = float(marg(p, THD).min())
    rho = rho_of(beta)
    L = cost(p)
    Lc = L*rho/(rho + mm) if mm < 0 else L    # exact repair by inflation
    return dict(p=p, L=float(Lc), Lraw=float(L), mm=mm)

def optimize_zig(beta, p0, fixr=None):
    def full(q):
        return np.r_[q[:8], fixr] if fixr is not None else q
    def cost(q):
        g = zig_geom(full(q))
        return g['L'] if g else 1e3
    def marg(q, thetas):
        g = zig_geom(full(q))
        if g is None:
            return -np.ones(len(thetas))
        return margins_zig(beta, g, thetas)
    n = 8 if fixr is not None else 9
    bounds = [(-2, 2)]*8 + ([] if fixr is not None else [(0.0, 0.9)])
    r = cutting_opt(beta, np.asarray(p0)[:n], cost,
                    lambda q: zig_geom(full(q)), marg, bounds)
    if zig_geom(full(r['p'])) is None or r['mm'] < -1e-5:
        return None
    r['p'] = full(r['p'])
    return r

def sweep_diagonal(betas):
    out = {}
    prev = None
    for bdeg in betas:
        beta = np.deg2rad(bdeg)
        Lc, ptsc = zee_classic(beta)
        # orientation scan of the classical zee
        rots = []
        for rot in np.arange(0, 180, 5):
            c, s = np.cos(np.deg2rad(rot)), np.sin(np.deg2rad(rot))
            R = np.array([[c, -s], [s, c]])
            g = dict(r=0.0, pts=ptsc @ R.T, L=Lc, arcs=[])
            rots.append((float(margins_zig(beta, g, THD).min()), rot))
        rots.sort(reverse=True)
        seeds = []
        for mmr, rot in rots[:2]:
            c, s = np.cos(np.deg2rad(rot)), np.sin(np.deg2rad(rot))
            R = np.array([[c, -s], [s, c]])
            zr = ptsc @ R.T
            seeds.append(np.r_[zr[0], zr[1], zr[2], zr[3], 0.0])
        if prev is not None:
            seeds.insert(0, prev)
        best = None
        for p0 in seeds:
            b0 = optimize_zig(beta, p0, fixr=0.0)
            if b0 and (best is None or b0['L'] < best['L']):
                best = b0
            for pp in ([b0['p']] if b0 else []) + [p0]:
                pr = pp.copy(); pr[8] = max(pr[8], 0.02)
                bf = optimize_zig(beta, pr, fixr=None)
                if bf and (best is None or bf['L'] < best['L']):
                    best = bf
        if best:
            prev = best['p'].copy()
            out[float(bdeg)] = dict(L=best['L'], mm=best['mm'],
                                    r=float(best['p'][8]),
                                    p=[float(v) for v in best['p']],
                                    Lclassic=float(Lc),
                                    zee_mm=float(rots[0][0]))
            print(f"{bdeg:6.2f}: Lzig={best['L']:.7f} classic={Lc:.7f} "
                  f"gain={Lc-best['L']:+.1e} r={best['p'][8]:.4f} "
                  f"mm={best['mm']:+.0e} zee_mm={rots[0][0]:+.1e}", flush=True)
        else:
            print(f"{bdeg:6.2f}: diagonal infeasible", flush=True)
    return out

# ---------------- boundary branches ---------------------------------------

def trap_margins(beta, q, thetas):
    a, b, c = q
    pts = np.array([[-b, 0], [b, 0], [a, c], [-a, c]])
    R0, R1, R2, sb, cb = rot_mats(beta)
    U = np.stack([np.cos(thetas), np.sin(thetas)], 1)
    def sup(M):
        return ((U @ M) @ pts.T).max(axis=1)
    return 2*cb*sup(R0) + sup(R1) + sup(R2) - 2*sb*cb

def trap_len(q):
    a, b, c = q
    return 2*b + 2*np.hypot(a - b, c)

def sweep_trapezoid(betas, seed_tbl):
    out = {}
    prev = None
    for bdeg in betas:
        beta = np.deg2rad(bdeg)
        u = np.tan(beta/2)
        seeds = [] if prev is None else [prev]
        if seed_tbl:
            k = min(seed_tbl, key=lambda row: abs(row[0] - u))
            seeds.append(np.array(k[1:]))
        best = None
        for p0 in seeds:
            r = cutting_opt(beta, p0, trap_len, None,
                            lambda q, th: trap_margins(beta, q, th),
                            [(0, 1.2), (0, 1.2), (0, 1.2)])
            if r['mm'] > -1e-5 and (best is None or r['L'] < best['L']):
                best = r
        if best:
            prev = best['p'].copy()
            out[float(bdeg)] = dict(L=best['L'], mm=best['mm'],
                                    p=[float(v) for v in best['p']])
            a, b, c = best['p']
            print(f"{bdeg:6.2f}: Ltrap={best['L']:.7f} a={a:.4f} b={b:.4f} "
                  f"c={c:.4f} mm={best['mm']:+.0e}", flush=True)
        else:
            print(f"{bdeg:6.2f}: trapezoid infeasible", flush=True)
    return out

def sweep_smooth(betas, p36):
    out = {}
    pz = np.asarray(p36, float)
    for bdeg in betas:
        beta = np.deg2rad(bdeg)
        def cost(q):
            return ZF.lengths(beta, q)[0]
        best = None
        for sc in ([1]*6, [1, 1.01, 1, 1.02, 0.95, 1.0],
                   [0.9, 1, 1, 1.05, 1, 0.98]):
            r = cutting_opt(beta, pz*np.array(sc), cost, None,
                            lambda q, th: ZF.margins(beta, q, th),
                            [(0.0, 1.0), (-np.pi/2, np.pi/2), (0, np.pi),
                             (0, np.pi), (0, 1.0), (0, 1.0)])
            if r['mm'] > -1e-5 and (best is None or r['L'] < best['L']):
                best = r
        if best is None:
            print(f"{bdeg:6.2f}: smooth infeasible", flush=True)
            continue
        px = best['p']
        # diameter-collapse detection: bar-only degenerate solution
        collapsed = (px[4] + px[5] < 5e-3 and
                     abs(2*px[0] - 2*np.cos(beta)) < 5e-2)
        pz = px.copy() if not collapsed else pz
        out[float(bdeg)] = dict(L=best['L'], mm=best['mm'],
                                p=[float(v) for v in px],
                                collapsed=bool(collapsed))
        print(f"{bdeg:6.2f}: Lsmooth={best['L']:.7f} q={px[0]:.4f} "
              f"gs={np.rad2deg(px[3]):5.2f} t+d={px[4]+px[5]:.4f} "
              f"mm={best['mm']:+.0e}{' [DIAM]' if collapsed else ''}", flush=True)
    return out

if __name__ == '__main__':
    # trapezoid seeds from the interval-theorem construction
    tbl = []
    try:
        BC = json.load(open('results/branch_curves.json'))
        for key in ('branch1', 'branch2'):
            for row in BC[key][::4]:
                tbl.append((row['u'], float(row['a']), float(row['b']),
                            float(row['c'])))
    except Exception as e:
        print("no branch_curves seeds:", e)

    print("== diagonal (zigzag) branch, r free ==")
    diag = sweep_diagonal(np.arange(60.0, 35.9, -1.0))

    print("\n== polygonal trapezoid branch ==")
    trap = sweep_trapezoid(np.arange(36.0, 46.1, 0.5), tbl)

    print("\n== smoothed (line-arc) boundary branch ==")
    p36 = [0.180832, np.deg2rad(-13.6894), 0.0, np.deg2rad(9.5234),
           0.071783, 0.291024]
    smooth = sweep_smooth(np.arange(36.0, 46.1, 0.5), p36)

    json.dump(dict(diag=diag, trap=trap, smooth=smooth),
              open('results/zigzag_branch.json', 'w'), indent=1)
    print("DONE")
