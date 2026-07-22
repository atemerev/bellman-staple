"""Exact characterization of the critical smoothed staple at the golden gnomon.

Contact structure (established numerically, margins sweep):
  (A) tangency ALONG the base-copy arc-image window  =>  the W-arc there has
      radius phi*r = rho (forcing r = sin(beta)) and center at the origin:
          phi*R0*O + R1*Va + R2*Vb = 0            (2 scalar equations)
  (B) isolated tangency at u = (-1,0):  h_W(u180) = rho   (1 equation)
  (C) KKT stationarity of the length within the 5-parameter template
      x = (ox, oy, psi1, psi2, d):  grad L = J^T lambda   (5 equations)
Newton at 60 digits, then PSLQ for minimal polynomials of the solved
quantities (tan-half-angles, coordinates, lengths).
"""
import json
import numpy as np
import mpmath as mp

mp.mp.dps = 70

S36 = mp.sin(mp.pi/5)      # sin 36
C36 = mp.cos(mp.pi/5)
PHI = 2*C36
RHO = 2*S36*C36            # sin 72
R = S36                    # curvature law: r = sin(beta), forced by (A)

R0 = mp.matrix([[0, -1], [1, 0]])
R1 = mp.matrix([[S36, C36], [-C36, S36]])
R2 = mp.matrix([[-S36, C36], [-C36, -S36]])
ROTS = [R0, R1, R2]
WTS = [PHI, mp.mpf(1), mp.mpf(1)]

def geom(x):
    ox, oy, p1, p2, d = x
    O = mp.matrix([ox, oy])
    T1 = mp.matrix([ox + R*mp.cos(p1), oy + R*mp.sin(p1)])
    t = T1[1] / mp.cos(p1)
    K = mp.matrix([T1[0] + t*mp.sin(p1), 0])
    T2 = mp.matrix([ox + R*mp.cos(p2), oy + R*mp.sin(p2)])
    E = mp.matrix([T2[0] - d*mp.sin(p2), T2[1] + d*mp.cos(p2)])
    return O, T1, T2, K, E, t

def verts_all(x):
    O, T1, T2, K, E, t = geom(x)
    Km = mp.matrix([-K[0], K[1]]); Em = mp.matrix([-E[0], E[1]])
    T1m = mp.matrix([-T1[0], T1[1]]); T2m = mp.matrix([-T2[0], T2[1]])
    return dict(K=K, E=E, Km=Km, Em=Em, T1=T1, T2=T2, T1m=T1m, T2m=T2m, O=O,
                Om=mp.matrix([-O[0], O[1]]))

def length(x):
    O, T1, T2, K, E, t = geom(x)
    return 2*(K[0] + t + R*(x[3]-x[2]) + x[4])

# ---- feature assignment at the two contacts (from the float solution) ----
xf = json.load(open('results/arc_template_opt.json'))['p']
xf = [mp.mpf(v) for v in [xf[0], xf[1], xf[3], xf[4], xf[5]]]

def copy_support_features(x, u, i):
    """(value, name) of the argmax feature of copy i at direction u."""
    V = verts_all(x)
    Rm, w = ROTS[i], WTS[i]
    best = (-mp.inf, None)
    for name in ('K', 'E', 'Km', 'Em', 'T1', 'T2', 'T1m', 'T2m'):
        p = w*(Rm*V[name])
        val = p[0]*u[0] + p[1]*u[1]
        if val > best[0]:
            best = (val, 'v:'+name)
    # arc caps (right arc normals [p1,p2]; left arc [pi-p2, pi-p1])
    ox, oy, p1, p2, d = x
    ru = mp.matrix([Rm[0,0]*u[0] + Rm[1,0]*u[1], Rm[0,1]*u[0] + Rm[1,1]*u[1]])
    ang = mp.atan2(ru[1], ru[0])
    for nm, ctr, lo, hi in (('arcR', V['O'], p1, p2),
                            ('arcL', V['Om'], mp.pi-p2, mp.pi-p1)):
        if lo <= ang <= hi:
            c = w*(Rm*ctr)
            val = c[0]*u[0] + c[1]*u[1] + w*R
            if val > best[0]:
                best = (val, 'a:'+nm)
    return best

u_win = mp.matrix([mp.cos(mp.mpf(81)*mp.pi/180), mp.sin(mp.mpf(81)*mp.pi/180)])
u180 = mp.matrix([-1, 0])
featsA = [copy_support_features(xf, u_win, i)[1] for i in range(3)]
featsB = [copy_support_features(xf, u180, i)[1] for i in range(3)]
print("window features:", featsA, "  u180 features:", featsB)

def support_of(x, i, feat, u):
    V = verts_all(x); Rm, w = ROTS[i], WTS[i]
    kind, name = feat.split(':')
    if kind == 'v':
        p = w*(Rm*V[name])
        return p[0]*u[0] + p[1]*u[1]
    ctr = V['O'] if name == 'arcR' else V['Om']
    c = w*(Rm*ctr)
    return c[0]*u[0] + c[1]*u[1] + w*R

def constraints(x):
    V = verts_all(x)
    # (A) center of the active W-arc at origin: phi*R0*O + sum of the two
    # vertex features on the window
    ctr = WTS[0]*(ROTS[0]*V['O'])
    for i in (1, 2):
        kind, name = featsA[i].split(':')
        assert kind == 'v'
        p = WTS[i]*(ROTS[i]*V[name])
        ctr = ctr + p
    # (B) support equality at u180
    hB = sum(support_of(x, i, featsB[i], u180) for i in range(3)) - RHO
    return [ctr[0], ctr[1], hB]

def grad(f, x, h=mp.mpf(10)**-25):
    g = []
    for i in range(5):
        xp = list(x); xm = list(x)
        xp[i] += h; xm[i] -= h
        g.append((f(xp) - f(xm)) / (2*h))
    return g

def F(z):
    x, lam = list(z[:5]), list(z[5:])
    C = constraints(x)
    gL = grad(length, x)
    Jc = [grad(lambda y, k=k: constraints(y)[k], x) for k in range(3)]
    stat = [gL[j] - sum(lam[k]*Jc[k][j] for k in range(3)) for j in range(5)]
    return C + stat

# multiplier seed by least squares
C0 = constraints(xf)
gL0 = grad(length, xf)
J0 = mp.matrix([grad(lambda y, k=k: constraints(y)[k], xf) for k in range(3)])
lam0 = mp.lu_solve(J0*J0.T, J0*mp.matrix(gL0))
z = list(xf) + [lam0[0], lam0[1], lam0[2]]
print("initial residual:", max(abs(v) for v in F(z)))

for it in range(80):
    Fz = F(z)
    J = mp.matrix(8, 8)
    h = mp.mpf(10)**-25
    for j in range(8):
        zp = list(z); zm = list(z)
        zp[j] += h; zm[j] -= h
        Fp, Fm = F(zp), F(zm)
        for i in range(8):
            J[i, j] = (Fp[i] - Fm[i]) / (2*h)
    dz = mp.lu_solve(J, mp.matrix([-v for v in Fz]))
    z = [z[i] + dz[i] for i in range(8)]
    if max(abs(v) for v in dz) < mp.mpf(10)**-62:
        break
print("newton done, residual:", max(abs(v) for v in F(z)), "iters", it)

x = z[:5]
ox, oy, p1, p2, d = x
O, T1, T2, K, E, t = geom(x)
L = length(x)
print("\nL  =", mp.nstr(L, 50))
print("q  =", mp.nstr(K[0], 45))
print("t  =", mp.nstr(t, 45))
print("d  =", mp.nstr(d, 45))
print("psi1 =", mp.nstr(p1*180/mp.pi, 30), "deg   psi2 =", mp.nstr(p2*180/mp.pi, 30))
print("O  =", mp.nstr(ox, 45), mp.nstr(oy, 45))
print("multipliers:", [mp.nstr(v, 12) for v in z[5:]])

# ---------------- PSLQ hunt -----------------------------------------------
targets = {
    'tan(psi1/2)': mp.tan(p1/2), 'tan(psi2/2)': mp.tan(p2/2),
    'sin(psi1)': mp.sin(p1), 'sin(psi2)': mp.sin(p2),
    'ox': ox, 'oy': oy, 'd': d, 'q': K[0], 't': t,
    'cos(dpsi)': mp.cos(p2-p1), 'E_x': E[0], 'E_y': E[1],
    'L_alg (L - 2R dpsi)': L - 2*R*(p2-p1),
}
for name, val in targets.items():
    found = None
    for deg in (1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24):
        pol = mp.findpoly(val, deg, maxcoeff=10**12, tol=mp.mpf(10)**-55)
        if pol:
            found = (deg, pol)
            break
    print(f"{name:22s} = {mp.nstr(val, 30)}  ->",
          f"deg {found[0]}: {found[1]}" if found else "no minpoly <= deg 24")

json.dump(dict(x=[mp.nstr(v, 60) for v in x],
               L=mp.nstr(L, 60), lam=[mp.nstr(v, 30) for v in z[5:]]),
          open('results/arc_exact_solution.json', 'w'), indent=1)
