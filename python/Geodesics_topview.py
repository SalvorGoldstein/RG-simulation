"""Top view (x-y plane) of N geodesics, Schwarzschild or Kerr.

Reads the files Geodesics<n>plus.dat (counter-clockwise, blue) and
Geodesics<n>minus.dat (clockwise, orange), columns: t r theta phi,
and the first line of dataforpy.txt (r_s).
"""
import glob
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# go to the project folder (this script is in a subfolder, hence '..')
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

SPIN = 0.0          # 0 = Schwarzschild, between 0 and 1 for Kerr (same value as in metric.h)
SAVE = None         # for example "topview.png" to also save the picture, None to only show

# ---------- data ----------
r_s = np.loadtxt('dataforpy.txt')[0]          # first line = r_s


def load(pattern):
    files = glob.glob(pattern)
    # sort by the number in the name (works for 3, 10, or 1.500000)
    files.sort(key=lambda f: float(re.search(r'Geodesics([\d.]+)', f).group(1)))
    paths = []
    for name in files:
        A = np.loadtxt(name, ndmin=2)
        if A.size == 0:
            continue
        r, th, ph = A[:, 1], A[:, 2], A[:, 3]
        paths.append((r * np.sin(th) * np.cos(ph),      # projection on the x-y plane
                      r * np.sin(th) * np.sin(ph)))
    return paths


ccw = load('data/Geodesics*plus.dat')      # counter-clockwise (blue)
cw = load('data/Geodesics*minus.dat')      # clockwise (orange)

if not ccw and not cw:
    raise SystemExit("No Geodesics*plus.dat or Geodesics*minus.dat file found. Run ./main first.")

# ---------- figure ----------
fig, ax = plt.subplots(figsize=(8, 8))


def draw(paths, cmap):
    cols = cmap(np.linspace(0.4, 0.95, max(len(paths), 1)))
    for (x, y), col in zip(paths, cols):
        ax.plot(x, y, color=col, lw=1.2)
        ax.plot(x[0], y[0], 'o', color=col, ms=4)       # starting point


draw(ccw, plt.cm.Blues)
draw(cw, plt.cm.Oranges)

# horizon (and ergosphere in the equatorial plane for Kerr)
phi = np.linspace(0, 2 * np.pi, 400)
r_h = r_s * (1 + np.sqrt(1 - SPIN**2)) / 2           # equals r_s for SPIN = 0
ax.fill(r_h * np.cos(phi), r_h * np.sin(phi), color='black', zorder=3)
handles = [Line2D([], [], color='black', marker='o', ls='', label='horizon')]
if SPIN > 0:
    ax.plot(r_s * np.cos(phi), r_s * np.sin(phi), 'r--', lw=1)
    handles.append(Line2D([], [], color='red', ls='--', label='ergosphere (equator)'))
if ccw:
    handles.append(Line2D([], [], color=plt.cm.Blues(0.7), label=f'counter-clockwise ({len(ccw)})'))
if cw:
    handles.append(Line2D([], [], color=plt.cm.Oranges(0.7), label=f'clockwise ({len(cw)})'))

ax.set_aspect('equal')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_title(f'{len(ccw) + len(cw)} geodesics, top view'
             + (f' (Kerr, spin = {SPIN})' if SPIN > 0 else ' (Schwarzschild)'))
ax.legend(handles=handles, loc='upper right')
ax.grid(alpha=0.3)

if SAVE:
    fig.savefig(SAVE, dpi=200, bbox_inches='tight')
plt.show()