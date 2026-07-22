"""Representative escape paths, one per regime, at critical placement
inside T(beta): zee@50, line-arc staple@41, Zalgalloid@30, caliper@26.
Emits TikZ coordinate blocks for the paper's gallery figure."""
import json
import warnings
import numpy as np
warnings.filterwarnings('ignore')
np.seterr(all='ignore')
import zalgalloid_family as ZF
from zigzag_branch import zee_classic

def tri_data(beta):
    sb, cb = np.sin(beta), np.cos(beta)
    N = np.array([[0., -1.], [sb, cb], [-sb, cb]])
    dvec = np.array([0., sb*cb, sb*cb])
    e = np.array([2*cb, 1., 1.])
    rho = 2*sb*cb
    return N, dvec, e, rho

def place(beta, pts, theta=None):
    """critical placement of path points inside T(beta); theta=None -> scan."""
    N, dvec, e, rho = tri_data(beta)
    def slack(th):
        c, s = np.cos(th), np.sin(th)
        R = np.array([[c, -s], [s, c]])
        Q = pts @ R.T
        h = np.array([ (Q @ n).max() for n in N ])
        return rho - e @ h, Q, h
    if theta is None:
        best = None
        for th in np.arange(0, 2*np.pi, np.deg2rad(0.02)):
            g = slack(th)[0]
            if best is None or g > best[0]:
                best = (g, th)
        theta = best[1]
    g, Q, h = slack(theta)
    t, *_ = np.linalg.lstsq(N, dvec - h, rcond=None)
    out = Q + t
    viol = max((out @ n).max() - dv for n, dv in zip(N, dvec))
    return out, float(g), float(viol)

def branch_path(beta, p, narc=13):
    """polyline of the boundary-type member (template frame, bar on y=0)."""
    q, nu0, gf, gs, tt, dd = p
    r = np.sin(beta)
    Kp, Cf, P1, P2, Cs, P3, E = (None,)*7
    rr, nu1, nu2, Kp, Cf, P1, P2, Cs, P3, E = ZF.hull_data(beta, p)
    right = [Kp]
    if gf > 1e-5:
        for a in np.linspace(nu0, nu1, narc)[1:]:
            right.append(Cf + r*np.array([np.cos(a), np.sin(a)]))
    if tt > 1e-9:
        right.append(P2)
    if gs > 1e-5:
        for a in np.linspace(nu1, nu2, narc)[1:]:
            right.append(Cs + r*np.array([np.cos(a), np.sin(a)]))
    right.append(E)
    right = np.array(right)
    left = right[::-1] * np.array([-1, 1])
    return np.vstack([left, right])

def fmt(pts, nd=4):
    return ' '.join(f'({x:.{nd}f},{y:.{nd}f})' for x, y in pts)

panels = []
zs = {round(r['base'], 2): r for r in json.load(open('results/zalgalloid_sweep.json'))}
zb = json.load(open('results/zigzag_branch.json'))

# --- zee @ 50 ---
beta = np.deg2rad(50)
L, pts = zee_classic(beta)
out, g, viol = place(beta, pts)
print(f"zee@50: L={L:.6f} slack={g:+.1e} viol={viol:+.1e}")
panels.append(dict(name='zee', beta=50.0, L=L, pts=fmt(out)))

# --- line-arc staple @ 41 ---
beta = np.deg2rad(41)
p = np.array(zb['smooth']['41.0']['p'])
pp = branch_path(beta, p)
out, g, viol = place(beta, -pp, theta=0.0)   # theta=180 pose == negate pts
Lb = zb['smooth']['41.0']['L']
print(f"staple@41: L={Lb:.6f} slack={g:+.1e} viol={viol:+.1e}")
panels.append(dict(name='staple', beta=41.0, L=Lb, pts=fmt(out)))

# --- Zalgalloid @ 30 ---
beta = np.deg2rad(30)
p = np.array(zs[30.0]['p'])
pp = branch_path(beta, p)
out, g, viol = place(beta, -pp, theta=0.0)
Lz = zs[30.0]['L']
print(f"zalgalloid@30: L={Lz:.6f} slack={g:+.1e} viol={viol:+.1e}")
panels.append(dict(name='zalgalloid', beta=30.0, L=Lz, pts=fmt(out)))

# --- caliper @ 26 (branch member == scaled caliper) ---
beta = np.deg2rad(26)
p = np.array(zs[26.0]['p'])
pp = branch_path(beta, p)
out, g, viol = place(beta, -pp, theta=0.0)
Lc = zs[26.0]['L']
print(f"caliper@26: L={Lc:.6f} slack={g:+.1e} viol={viol:+.1e}")
panels.append(dict(name='caliper', beta=26.0, L=Lc, pts=fmt(out)))

json.dump(panels, open('results/gallery_panels.json', 'w'), indent=1)
for pn in panels:
    b = np.deg2rad(pn['beta'])
    print(f"\n%% {pn['name']} beta={pn['beta']} L={pn['L']:.6f}")
    print(f"%% triangle: (-{np.cos(b):.5f},0) -- ({np.cos(b):.5f},0) -- (0,{np.sin(b):.5f})")
    print(pn['pts'])
