"""Zalgalloid continuation figure: six optimized members, apex-angle
labels, exact-anchor tag on the 108-degree panel, dashed caliper overlay
on the last panel.  Sized 1:1 for \textwidth inclusion."""
import json
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

C_BR = '#eb6834'
C_TRI = '0.40'
C_CAL = '0.15'
C_ANCHOR = '#4a3aa7'

panels = json.load(open('results/zalgalloid_panels_tikz.json'))
BETAS = [36.0, 33.0, 31.0, 29.0, 28.0, 27.25]
LVALS = [1.282676, 1.218281, 1.168117, 1.104089, 1.069529, 1.043170]

# caliper overlay for the last panel, from the archived TikZ figure
src = open('paper/staple_v3_backup.tex').read()
i = src.find('teal!75!black,densely dashed] plot coordinates {')
j = src.find('};', i)
cal = np.array(re.findall(r'\(([-\d.]+),([-\d.]+)\)', src[i:j]), float)

plt.rcParams.update({'font.size': 9, 'mathtext.fontset': 'cm'})
fig, axes = plt.subplots(2, 3, figsize=(6.3, 2.62), dpi=300)
fig.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.02,
                    wspace=0.04, hspace=0.0)

for k, (bdeg, L, ax) in enumerate(zip(BETAS, LVALS, axes.ravel())):
    beta = np.deg2rad(bdeg)
    cb, sb = np.cos(beta), np.sin(beta)
    ax.plot([-cb, cb, 0, -cb], [0, 0, sb, 0], color=C_TRI, lw=1.0)
    pts = np.array(panels[str(bdeg)])
    ax.plot(pts[:, 0], pts[:, 1], color=C_BR, lw=1.9,
            solid_joinstyle='round', solid_capstyle='round')
    if k == 5:
        ax.plot(cal[:, 0], cal[:, 1], color=C_CAL, lw=1.1,
                ls=(0, (2.6, 2.0)))
    alpha = 180 - 2*bdeg
    line1 = f'$\\alpha={alpha:g}^\\circ$'
    line2 = f'$L={L:.6f}$'
    if k == 5:
        line2 = f'$L={L:.6f}$ $=$ caliper (dashed)'
    ax.text(0, -0.115, line1, ha='center', va='top', fontsize=8.6,
            transform=ax.transData)
    ax.text(0, -0.255, line2, ha='center', va='top', fontsize=7.4,
            color='0.30')
    if k == 0:
        ax.plot([-0.87], [0.56], marker='*', ms=8, color=C_ANCHOR)
        ax.text(-0.79, 0.56, 'exact anchor', fontsize=7.0, color=C_ANCHOR,
                ha='left', va='center')
    ax.set_xlim(-0.95, 0.95)
    ax.set_ylim(-0.36, 0.68)
    ax.set_aspect('equal')
    ax.axis('off')

fig.savefig('results/fig_zalgalloid_morph.png')
print('wrote results/fig_zalgalloid_morph.png')
