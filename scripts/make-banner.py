#!/usr/bin/env python3
# make-banner.py — generate the GitHub social-preview / README banner.
# Usage (from the repo root, any python3 with matplotlib):
#   python scripts\make-banner.py
# Output: assets\banner.png (1280x640)
import os
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

ASSET = os.path.join(os.path.dirname(__file__), '..', 'assets', 'banner.png')
os.makedirs(os.path.dirname(ASSET), exist_ok=True)

# Chinese font (optional; fall back to English-only tagline)
zh_family = None
for cand in ('C:/Windows/Fonts/msyh.ttc', 'C:/Windows/Fonts/msyh.ttf',
             'C:/Windows/Fonts/simhei.ttf'):
    if os.path.exists(cand):
        try:
            font_manager.fontManager.addfont(cand)
            zh_family = font_manager.FontProperties(fname=cand).get_name()
            break
        except Exception:
            pass

W, H = 12.8, 6.4
fig, ax = plt.subplots(figsize=(W, H), dpi=100)
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis('off')

# gradient background: deep navy -> teal-ish blue
grad = np.linspace(0, 1, 256).reshape(1, -1)
top = np.array([13, 22, 38]) / 255.0      # #0D1626
bot = np.array([24, 58, 96]) / 255.0      # #183A60
bg = np.zeros((256, 256, 3))
for i in range(3):
    bg[:, :, i] = top[i] * (1 - grad) + bot[i] * grad
ax.imshow(bg, extent=[0, W, 0, H], aspect='auto', origin='lower', zorder=0)

# subtle decorative "molecule": nodes + bonds in accent teal
accent = '#3EC6A8'
nodes = [(10.3, 3.4), (10.85, 4.05), (11.45, 3.55), (12.1, 4.2)]
bonds = [(0, 1), (1, 2), (2, 3), (0, 2)]
for a, b in bonds:
    ax.plot([nodes[a][0], nodes[b][0]], [nodes[a][1], nodes[b][1]],
            color=accent, lw=2.2, alpha=0.75, zorder=1)
for x, y in nodes:
    ax.scatter(x, y, s=340, color=accent, alpha=0.9, zorder=2, edgecolors='none')

# title
ax.text(0.9, 4.35, 'dsh-bioinfo', fontsize=64, fontweight='bold', color='white',
        va='center', ha='left', zorder=3)
ax.text(0.95, 3.28, '生信模式', fontsize=34, color='#D8E4F2', va='center', ha='left',
        zorder=3, fontfamily=(zh_family or 'sans-serif'))
ax.text(0.95, 2.45, 'Bioinformatics agent preset for DeepSeek Harness',
        fontsize=21, color='#D8E4F2', va='center', ha='left', zorder=3)
ax.text(0.95, 1.55,
        'AF2-Multimer  ·  ESMFold  ·  AutoDock Vina  ·  MM-GBSA / MD  ·  '
        'TM-score QA  ·  Virtual Screening',
        fontsize=15, color='#9FB6CC', va='center', ha='left', zorder=3)

fig.savefig(ASSET, dpi=100, facecolor=fig.get_facecolor())
print('written', os.path.abspath(ASSET))
