"""The candidate upper-bound frontier from the equilateral triangle to the
flat isosceles limit: apex angle 60 <= alpha < 180 degrees.

Design: main panel shows the envelope plus, in thin dashes, everything that
makes the two events legible — the zee and branch continuations past the
crossing (transversality) and the scaled caliper approaching the branch
(tangential merge).  Events and rigorous status live in a separate lane so
no text competes with the curves.  Sized 1:1 for \textwidth inclusion.
"""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ZETA = 2.2782916414368529
BX = 42.2867          # zee/branch crossing (transversal)
BSTAR = 27.62         # tangential junction with the caliper
AX = 180.0 - 2.0*BX
ASTAR = 180.0 - 2.0*BSTAR

C_ZEE = '#2a78d6'     # series 1 (validated pair)
C_BR = '#eb6834'      # series 2
C_ZEE_INK = '#174d91' # label ink, darker steps of the same hues
C_BR_INK = '#9c3d14'
C_CAL = '0.45'
C_CROSS = '#c1272d'
C_ANCHOR = '#4a3aa7'

zs = json.load(open('results/zalgalloid_sweep.json'))
zb = json.load(open('results/zigzag_branch.json'))
cf = json.load(open('results/crossing_fine.json'))

branch = {}
for r in zs:
    branch[round(r['base'], 2)] = r['L']
for k, v in zb['smooth'].items():
    if not v['collapsed'] and float(k) <= 45.0:
        branch[round(float(k), 2)] = v['L']
for b, lz, ls, gs in cf['fine']:
    branch[round(b, 2)] = ls
for b, gs, ls in cf['gs_collapse']:
    branch[round(b, 2)] = ls
branch[45.0] = float(np.sqrt(2))
bB = np.array(sorted(branch))
LB = np.array([branch[b] for b in bB])
aB = 180.0 - 2.0*bB
order = np.argsort(aB)
aB, LB = aB[order], LB[order]

def zee_from_beta(beta):
    sb, cb = np.sin(beta), np.cos(beta)
    return 6*sb*cb/np.sqrt(1 + 8*sb**2)

def beta_from_apex(alpha_deg):
    return np.deg2rad((180.0 - alpha_deg)/2.0)

afull = np.linspace(60.0, 180.0, 900)
caliper = ZETA*np.cos(np.deg2rad(afull)/2.0)

az = np.linspace(60.0, 108.0, 400)
Lz = zee_from_beta(beta_from_apex(az))
Lx = 1.3891993
Ls = ZETA*np.sin(np.deg2rad(BSTAR))

plt.rcParams.update({'font.size': 9, 'mathtext.fontset': 'cm',
                     'axes.linewidth': 0.7})
fig = plt.figure(figsize=(6.3, 4.75), dpi=300)
gsp = fig.add_gridspec(2, 1, height_ratios=[3.55, 1.0], hspace=0.045)
ax = fig.add_subplot(gsp[0])
band = fig.add_subplot(gsp[1], sharex=ax)

ax.axvspan(60, AX, color=C_ZEE, alpha=0.045)
ax.axvspan(AX, 180, color=C_BR, alpha=0.05)
for xevent in (AX, 108.0, ASTAR):
    ax.axvline(xevent, color='0.62', lw=0.6, ls=(0, (1.5, 2.4)), zorder=0)

# --- everything the events need, in thin dashes -------------------------
zc = az >= AX                                     # zee continuation
ax.plot(az[zc], Lz[zc], color=C_ZEE, lw=1.0, ls=(0, (4, 2.6)), alpha=0.9)
bc = aB <= AX                                     # branch continuation
ax.plot(aB[bc], LB[bc], color=C_BR, lw=1.0, ls=(0, (4, 2.6)), alpha=0.9)
mcalref = (caliper <= 1.505) & (afull <= ASTAR)   # caliper approaching
ax.plot(afull[mcalref], caliper[mcalref], color=C_CAL, lw=1.1,
        ls=(0, (5, 2.2)))

# --- the frontier itself -------------------------------------------------
zf = az <= AX
ax.plot(az[zf], Lz[zf], color=C_ZEE, lw=2.6, solid_capstyle='round')
bf = (aB > AX) & (aB < ASTAR)
abranch = np.concatenate(([AX], aB[bf], [ASTAR]))
Lbranch = np.concatenate(([Lx], LB[bf], [Ls]))
ax.plot(abranch, Lbranch, color=C_BR, lw=2.6, solid_capstyle='round')
mcal = afull >= ASTAR
ax.plot(afull[mcal], caliper[mcal], color=C_BR, lw=2.6,
        solid_capstyle='round')

ax.plot([AX], [Lx], marker='D', ms=6, color=C_CROSS, zorder=6)
ax.plot([108], [1.2826799], marker='*', ms=9, color=C_ANCHOR, zorder=6)
ax.plot([ASTAR], [Ls], marker='o', ms=6.5, mfc='white', mec=C_BR, mew=1.5,
        zorder=6)

ax.text(73, zee_from_beta(beta_from_apex(73)) + 0.055, 'zee',
        fontsize=10, color=C_ZEE_INK, ha='center')
ax.text(149.5, ZETA*np.cos(np.deg2rad(149.5)/2.0) + 0.075,
        'Zalgalloid\u2013Zalgaller branch', fontsize=10, color=C_BR_INK,
        ha='center', rotation=-31, rotation_mode='anchor')
ax.text(114.0, 1.318, r'scaled caliper $\zeta\sin\beta$', fontsize=8,
        color='0.32', ha='left', rotation=-30, rotation_mode='anchor')

ax.set_xlim(60, 180)
ax.set_ylim(0, 1.5)
ax.set_yticks(np.arange(0, 1.51, 0.25))
ax.set_ylabel('escape-path length (legs $= 1$)', fontsize=9)
ax.tick_params(axis='x', bottom=False, labelbottom=False)
ax.grid(axis='y', alpha=0.22, lw=0.5)
secax = ax.secondary_xaxis(
    'top', functions=(lambda alpha: (180.0-alpha)/2.0,
                      lambda beta: 180.0-2.0*beta))
secax.set_xticks([60, 45, 36, 30, 15])
secax.set_xlabel(r'base angle $\beta=(180^\circ-\alpha)/2$', fontsize=9)

# --- event + status lane -------------------------------------------------
band.set_ylim(0, 1)
band.plot([AX], [0.80], marker='D', ms=5, color=C_CROSS, clip_on=False)
band.text(AX - 2.6, 0.80, r'$95.43^\circ$ zee/branch crossing',
          fontsize=7.8, ha='right', va='center')
band.plot([108], [0.52], marker='*', ms=7, color=C_ANCHOR, clip_on=False)
band.text(110.0, 0.52, '$108^\\circ$ exact line\u2013arc anchor',
          fontsize=7.8, ha='left', va='center')
band.plot([ASTAR], [0.80], marker='o', ms=5.5, mfc='white', mec=C_BR,
          mew=1.3, clip_on=False)
band.text(ASTAR + 2.6, 0.80, r'$124.8^\circ$ smooth merge into the caliper',
          fontsize=7.8, ha='left', va='center')

acert = [180 - 2*42.1633, 180 - 2*31.0719]
band.plot([60, 90], [0.10, 0.10], color=C_ZEE, lw=3.4,
          solid_capstyle='butt')
band.text(75, 0.185, 'zee proved optimal', fontsize=7.8, ha='center',
          color=C_ZEE_INK, va='bottom')
band.plot(acert, [0.10, 0.10], color=C_BR, lw=3.4, solid_capstyle='butt')
band.text(np.mean(acert), 0.185, 'exact staple upper bounds',
          fontsize=7.8, ha='center', color=C_BR_INK, va='bottom')

band.set_xticks(np.arange(60, 181, 15))
band.set_xlabel(r'apex angle $\alpha$ (degrees), equal legs $=1$')
band.set_yticks([])
for side in ('left', 'right', 'top'):
    band.spines[side].set_visible(False)
band.tick_params(axis='x', pad=2.5, labelsize=8.5)
fig.subplots_adjust(left=0.10, right=0.972, top=0.905, bottom=0.105)
fig.savefig('results/fig_regimes.png')
print('wrote results/fig_regimes.png')

# thinned branch coordinates for pgfplots
sel = []
for b in bB:
    if b <= 42.0 and (abs(b*2 - round(b*2)) < 1e-9):
        sel.append(b)
    elif 42.0 < b:
        sel.append(b)
coords = ' '.join(f'({b:.2f},{branch[b]:.6f})' for b in sel)
open('results/regimes_pgf.txt', 'w').write(coords + '\n')
print('pgf coords:', len(sel), 'points -> results/regimes_pgf.txt')
json.dump(dict(branch={f'{b:.2f}': branch[b] for b in sorted(branch)},
               crossing=cf['crossing'], bstar=BSTAR, astar=ASTAR,
               arc_on=42.88, rect_endpoint=[45.0, float(np.sqrt(2))]),
          open('results/regime_diagram.json', 'w'), indent=1)
print('DONE')
