"""Assemble the isosceles regime diagram: data JSON, PNG figure, and the
pgfplots coordinate block for the paper."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ZETA = 2.2782916414368529
BX = 42.2867          # zee/staple crossing (first-order)
BSTAR = 27.62         # tangential junction with the caliper (second-order)

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

def zee(beta):
    sb, cb = np.sin(beta), np.cos(beta)
    return 6*sb*cb/np.sqrt(1 + 8*sb**2)

bz = np.linspace(np.deg2rad(36), np.deg2rad(60), 200)
bfull = np.linspace(np.deg2rad(26), np.deg2rad(60), 400)

fig, ax = plt.subplots(figsize=(9.5, 5.6), dpi=170)
ax.axvspan(26, BSTAR, color='teal', alpha=0.07)
ax.axvspan(BSTAR, BX, color='darkorange', alpha=0.08)
ax.axvspan(BX, 60, color='green', alpha=0.07)

ax.plot(np.rad2deg(bfull), 2*np.cos(bfull), color='0.55', lw=1.0, ls=(0, (2, 2)),
        label='diameter $2\\cos\\beta$')
ax.plot(np.rad2deg(bfull), ZETA*np.sin(bfull), color='teal', lw=1.0, alpha=0.65)
mcal = np.rad2deg(bfull) <= BSTAR
ax.plot(np.rad2deg(bfull)[mcal], (ZETA*np.sin(bfull))[mcal], color='teal',
        lw=2.6, label='scaled caliper $\\zeta\\sin\\beta$')
mzee = np.rad2deg(bz) >= BX
ax.plot(np.rad2deg(bz)[~mzee], zee(bz)[~mzee], color='green', lw=1.0,
        ls='--', alpha=0.8)
ax.plot(np.rad2deg(bz)[mzee], zee(bz)[mzee], color='green', lw=2.6,
        label='zee (proven optimal on $[45^\\circ,60^\\circ]$)')
mb = bB <= BX
ax.plot(bB[~mb], LB[~mb], color='darkorange', lw=1.0, ls='--', alpha=0.85)
ax.plot(bB[mb], LB[mb], color='darkorange', lw=2.6,
        label='line--arc branch (staple $\\to$ Zalgalloid)')
ax.plot([45], [np.sqrt(2)], marker='o', ms=6, mfc='white', mec='darkorange',
        zorder=6)
ax.annotate('exact rectangle endpoint\n$L=\\sqrt{2}$ = diameter tie',
            xy=(45, np.sqrt(2)), xytext=(46.6, 1.405), fontsize=8.5,
            arrowprops=dict(arrowstyle='-', lw=0.7, color='0.4'))

Lx = 1.3891993
ax.plot([BX], [Lx], marker='D', ms=6, color='crimson', zorder=6)
ax.annotate('first-order crossing $\\beta_\\times=42.287^\\circ$\n'
            'kink = hardest isosceles, $L=1.38920$\n'
            'slopes $-0.016$ / $+0.013$ per degree',
            xy=(BX, Lx), xytext=(46.8, 1.30), fontsize=8.5,
            arrowprops=dict(arrowstyle='->', lw=0.9, color='crimson'))
Ls = ZETA*np.sin(np.deg2rad(BSTAR))
ax.plot([BSTAR], [Ls], marker='o', ms=7, mfc='white', mec='teal', mew=1.6,
        zorder=6)
ax.annotate('tangential junction $\\beta^{*}\\!\\approx27.6^\\circ$\n'
            '(branch merges into the caliper)',
            xy=(BSTAR, Ls), xytext=(27.4, 1.155), fontsize=8.5,
            arrowprops=dict(arrowstyle='->', lw=0.9, color='teal'))
ax.plot([36], [1.2826760], marker='v', ms=6, color='purple', zorder=6)
ax.annotate('golden gnomon:\ncertified $\\leq 1.2826799$',
            xy=(36, 1.28268), xytext=(33.2, 1.335), fontsize=8.5,
            arrowprops=dict(arrowstyle='-', lw=0.7, color='purple'))
ax.plot([31.07, 42.16], [0.985, 0.985], color='darkorange', lw=3, alpha=0.9,
        solid_capstyle='butt')
ax.text(36.6, 0.992, 'Theorem 1: certified interval $[31.07^\\circ,42.16^\\circ]$',
        fontsize=8.5, ha='center', color='saddlebrown')

for x, lab, col in ((26.8, 'CALIPER', 'teal'), (34.8, 'STAPLE / LINE--ARC',
                    'darkorange'), (51.5, 'ZEE', 'green')):
    ax.text(x, 1.445, lab, fontsize=10, fontweight='bold', color=col,
            ha='center', va='top')

ax.set_xlim(26, 60); ax.set_ylim(0.965, 1.455)
ax.set_xlabel('base angle $\\beta$ (degrees), legs $=1$')
ax.set_ylabel('escape-path length')
ax.legend(loc='lower right', fontsize=8.5, framealpha=0.95)
ax.grid(alpha=0.25)
fig.tight_layout()
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
json.dump(dict(branch={f'{b:.2f}': branch[b] for b in bB},
               crossing=cf['crossing'], bstar=BSTAR,
               arc_on=42.88, rect_endpoint=[45.0, float(np.sqrt(2))]),
          open('results/regime_diagram.json', 'w'), indent=1)
print('DONE')
