"""Fine localization of (i) the zee/staple regime crossing, (ii) the arc
turn-off point inside the boundary branch, (iii) the branch endpoint at 45.
"""
import json
import warnings
import numpy as np
warnings.filterwarnings('ignore')
np.seterr(all='ignore')
import zigzag_branch as ZB
import zalgalloid_family as ZF

D = json.load(open('results/zigzag_branch.json'))

def zee_opt(bdeg):
    """optimized diagonal length (classical zee refined)."""
    beta = np.deg2rad(bdeg)
    keys = sorted(D['diag'], key=lambda k: abs(float(k) - bdeg))
    p0 = np.array(D['diag'][keys[0]]['p'])
    b = ZB.optimize_zig(beta, p0, fixr=0.0)
    return b['L'] if b else None

def smooth_opt(bdeg, seed):
    beta = np.deg2rad(bdeg)
    best = None
    for sc in ([1]*6, [1, 1.005, 1, 1.01, 0.99, 1]):
        r = ZB.cutting_opt(beta, np.array(seed)*np.array(sc),
                           lambda q: ZF.lengths(beta, q)[0], None,
                           lambda q, th: ZF.margins(beta, q, th),
                           [(0.0, 1.0), (-np.pi/2, np.pi/2), (0, np.pi),
                            (0, np.pi), (0, 1.0), (0, 1.0)])
        if r['mm'] > -1e-5 and (best is None or r['L'] < best['L']):
            best = r
    return best

rows = []
seed = D['smooth']['42.0']['p']
print("beta    L_zee       L_smooth    gs(deg)  diff")
for bdeg in np.arange(42.05, 42.56, 0.05):
    lz = zee_opt(bdeg)
    bs = smooth_opt(bdeg, seed)
    seed = bs['p']
    gs = np.rad2deg(bs['p'][3])
    rows.append((float(bdeg), float(lz), float(bs['L']), float(gs)))
    print(f"{bdeg:6.2f} {lz:.8f} {bs['L']:.8f} {gs:7.3f} {bs['L']-lz:+.6f}",
          flush=True)

# crossing by local linear fit of the difference
import numpy.polynomial.polynomial as P
arr = np.array(rows)
d = arr[:, 2] - arr[:, 1]
i = np.argmin(np.abs(d))
sel = slice(max(0, i-2), i+3)
cf = np.polyfit(arr[sel, 0], d[sel], 2)
roots = np.roots(cf)
bx = [r for r in roots if arr[0, 0] - 0.2 < r < arr[-1, 0] + 0.2][0]
sz = np.polyfit(arr[sel, 0], arr[sel, 1], 2)
ss = np.polyfit(arr[sel, 0], arr[sel, 2], 2)
szp = np.polyval(np.polyder(sz), bx)
ssp = np.polyval(np.polyder(ss), bx)
Lx = np.polyval(sz, bx)
print(f"\nCROSSING: beta_x = {bx:.4f} deg, L = {Lx:.7f}")
print(f"slopes per degree: zee {szp:+.5f}, staple {ssp:+.5f}, gap {ssp-szp:.5f}")

# arc turn-off: continue smoothed branch 42.5 -> 43.2 finely
print("\narc span collapse:")
gs_rows = []
for bdeg in np.arange(42.50, 43.21, 0.1):
    bs = smooth_opt(bdeg, seed)
    if bs is None:
        continue
    seed = bs['p']
    gs = np.rad2deg(bs['p'][3])
    gs_rows.append((float(bdeg), float(gs), float(bs['L'])))
    print(f"{bdeg:6.2f}: gs={gs:7.4f} deg  L={bs['L']:.8f}", flush=True)

# rectangle endpoint at 45 deg: exact-candidate check
beta = np.deg2rad(45)
b = 1 - np.sqrt(2)/2; c = np.sqrt(2) - 1
mm = float(ZB.trap_margins(beta, [b, b, c], ZB.THD).min())
print(f"\nrect endpoint check at 45deg: a=b=1-1/sqrt2, c=sqrt2-1 -> "
      f"L={ZB.trap_len([b,b,c]):.10f} (sqrt2={np.sqrt(2):.10f}) mm={mm:+.2e}")

json.dump(dict(fine=rows, crossing=dict(beta=float(bx), L=float(Lx),
                                        slope_zee=float(szp),
                                        slope_staple=float(ssp)),
               gs_collapse=gs_rows),
          open('results/crossing_fine.json', 'w'), indent=1)
print("DONE")
