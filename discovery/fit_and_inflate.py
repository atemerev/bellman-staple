"""Fit rational-coefficient polynomials to the three staple branch curves,
inflate by (1+eps), and float-validate the support-form escape certificate
plus all length comparisons on a dense grid.

Support form: W = sum_i w_i R_i H contains disk(rho) iff for each of the 12
candidate edge-normal directions nu (= R_i g, g an edge normal of H),
h_W(nu) >= rho*|nu|.  We lower-bound h_W(nu) by fixed per-copy maximizer
vertices; validity of the lower bound is automatic.
"""
import itertools, json
from fractions import Fraction
import numpy as np

DATA = json.load(open('results/branch_curves.json'))
EPS = Fraction(1, 50000)          # scaling inflation
DEN = 10**9                       # coefficient denominator

SEGS = [
    ('branch1', Fraction(278, 1000), Fraction(3530, 10000), 6),
    ('vertex',  Fraction(3530, 10000), Fraction(353399, 10**6), 2),
    ('branch2', Fraction(353399, 10**6), Fraction(3855, 10000), 6),
]

# per-claim u-intervals (comparisons hold on these; ESC = escape everywhere)
CLAIMS = {
    'ESC': (Fraction(278, 1000), Fraction(3855, 10000)),
    'SQA': (Fraction(278, 1000), Fraction(73, 200)),      # vs base-aligned square
    'SQB': (Fraction(9, 25),     Fraction(3855, 10000)),  # vs leg-aligned square
    'ZEE': (Fraction(278, 1000), Fraction(3855, 10000)),
    'CAL': (Fraction(141, 500),  Fraction(3855, 10000)),
    'DIA': (Fraction(278, 1000), Fraction(3855, 10000)),
}

def fit_segment(name, lo, hi, deg):
    # fit on ALL solved points of the branch (which extend slightly beyond the
    # segment) so that [lo, hi] is interior to the data — no extrapolation
    recs = DATA[name]
    us = np.array([r['u'] for r in recs])
    out = {}
    for key in 'abc':
        ys = np.array([float(r[key]) for r in recs])
        cf = np.polynomial.polynomial.polyfit(us, ys, deg)
        resid = np.abs(np.polynomial.polynomial.polyval(us, cf) - ys).max()
        # rationalize and inflate
        cf_r = [Fraction(round(v * DEN), DEN) * (1 + EPS) for v in cf]
        out[key] = cf_r
        out[key + '_resid'] = resid
    return out, len(recs)

# ---------------- float evaluation of certificates --------------------
def polyval_fr(cf, u):
    return float(sum(float(c) * u**k for k, c in enumerate(cf)))

def tri(u):
    den = 1 + u*u
    return 2*u/den, (1 - u*u)/den     # sin beta, cos beta

def support_margins(u, a, b, c):
    sb, cb = tri(u)
    H = np.array([[-b, 0], [b, 0], [a, c], [-a, c]])
    nws = [(np.array([0.0, -1.0]), 2*cb), (np.array([sb, cb]), 1.0),
           (np.array([-sb, cb]), 1.0)]
    Rs = [np.array([[n[0], n[1]], [-n[1], n[0]]]) for n, _ in nws]
    ws = [w for _, w in nws]
    rho = 2*sb*cb
    # H edge normals (unnormalized): base, top, right leg, left leg
    gs = [np.array([0.0, -1.0]), np.array([0.0, 1.0]),
          np.array([c, -(a - b)]), np.array([-c, -(a - b)])]
    margins = []
    for i, R in enumerate(Rs):
        for g in gs:
            nu = R @ g
            hW = sum(ws[m] * (H @ (Rs[m].T @ nu)).max() for m in range(3))
            margins.append(hW - rho*np.linalg.norm(nu))
    return np.array(margins)

def support_margins_fixed(u, a, b, c, choices):
    """lower-bound margins with FIXED per-(direction,copy) maximizer vertices,
    exactly as the symbolic certificate will use."""
    sb, cb = tri(u)
    H = np.array([[-b, 0], [b, 0], [a, c], [-a, c]])
    nws = [(np.array([0.0, -1.0]), 2*cb), (np.array([sb, cb]), 1.0),
           (np.array([-sb, cb]), 1.0)]
    Rs = [np.array([[n[0], n[1]], [-n[1], n[0]]]) for n, _ in nws]
    ws = [w for _, w in nws]
    rho = 2*sb*cb
    gs = [np.array([0.0, -1.0]), np.array([0.0, 1.0]),
          np.array([c, -(a - b)]), np.array([-c, -(a - b)])]
    margins = []
    j = 0
    for i, R in enumerate(Rs):
        for g in gs:
            nu = R @ g
            hW = sum(ws[m] * (H[choices[j][m]] @ (Rs[m].T @ nu)) for m in range(3))
            margins.append(hW - rho*np.linalg.norm(nu))
            j += 1
    return np.array(margins)

def maximizer_indices(u, a, b, c):
    """which H-vertex attains the support for each (copy m, direction j)."""
    sb, cb = tri(u)
    H = np.array([[-b, 0], [b, 0], [a, c], [-a, c]])
    nws = [(np.array([0.0, -1.0]), 2*cb), (np.array([sb, cb]), 1.0),
           (np.array([-sb, cb]), 1.0)]
    Rs = [np.array([[n[0], n[1]], [-n[1], n[0]]]) for n, _ in nws]
    gs = [np.array([0.0, -1.0]), np.array([0.0, 1.0]),
          np.array([c, -(a - b)]), np.array([-c, -(a - b)])]
    out = []
    for i, R in enumerate(Rs):
        for g in gs:
            nu = R @ g
            out.append(tuple(int(np.argmax(H @ (Rs[m].T @ nu))) for m in range(3)))
    return out

def L_staple(u, a, b, c):
    return 2*b + 2*np.hypot(a - b, c)

ZETA_LO = 2.2782916

def comparisons(u, a, b, c):
    sb, cb = tri(u)
    L = L_staple(u, a, b, c)
    sqA = 6*sb*cb/(sb + 2*cb)
    sqB = 3*sb/(sb + cb)
    zee = 6*sb*cb/np.sqrt(1 + 8*sb*sb)
    cal = ZETA_LO*sb
    dia = 2*cb
    return dict(SQA=sqA - L, SQB=sqB - L, ZEE=zee - L, CAL=cal - L, DIA=dia - L)

if __name__ == '__main__':
    construction = {}
    for name, lo, hi, deg in SEGS:
        fit, npts = fit_segment(name, lo, hi, deg)
        construction[name] = dict(
            lo=[lo.numerator, lo.denominator], hi=[hi.numerator, hi.denominator],
            a=[[f.numerator, f.denominator] for f in fit['a']],
            b=[[f.numerator, f.denominator] for f in fit['b']],
            c=[[f.numerator, f.denominator] for f in fit['c']])
        print(f"{name}: {npts} pts deg={deg} resid a={fit['a_resid']:.1e} "
              f"b={fit['b_resid']:.1e} c={fit['c_resid']:.1e}")

        # fixed maximizer choice at segment midpoint
        umid = (float(lo) + float(hi))/2
        amid = polyval_fr(fit['a'], umid); bmid = polyval_fr(fit['b'], umid)
        cmid = polyval_fr(fit['c'], umid)
        choices = maximizer_indices(umid, amid, bmid, cmid)
        construction[name]['maximizers'] = choices

        # dense float validation with the FIXED choices (mirrors exact proof),
        # per claim on (claim interval) ∩ (segment)
        for cname, (clo, chi) in CLAIMS.items():
            plo, phi = max(lo, clo), min(hi, chi)
            if plo >= phi:
                continue
            uu = np.linspace(float(plo), float(phi), 3000)
            worst, worst_u = np.inf, None
            for u in uu:
                a = polyval_fr(fit['a'], u); b = polyval_fr(fit['b'], u)
                c = polyval_fr(fit['c'], u)
                if cname == 'ESC':
                    m = support_margins_fixed(u, a, b, c, choices).min()
                else:
                    m = comparisons(u, a, b, c)[cname]
                if m < worst:
                    worst, worst_u = m, u
            flag = 'OK' if worst > 0 else 'FAIL'
            print(f"   {cname} on [{float(plo):.4f},{float(phi):.4f}]: "
                  f"min {worst:+.3e} at u={worst_u:.5f}  {flag}")
    json.dump(construction, open('results/interval_construction.json', 'w'), indent=1)
    print("saved results/interval_construction.json")
