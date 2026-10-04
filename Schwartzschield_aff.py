import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation

# ---------- quality preset: "draft", "medium", "high" ----------
QUALITY = "high"

PRESETS = {
#            frames fps  dpi  inches  points  sphere_n
"draft":  dict(frames=120, fps=20, dpi=70,  size=5.0,  points=1500, n=12),
"medium": dict(frames=240, fps=30, dpi=100, size=7.2,  points=3000, n=24),   # 720x720
"high":   dict(frames=360, fps=30, dpi=120, size=9.0,  points=5000, n=40),   # 1080x1080
}
P = PRESETS[QUALITY]
N_FRAMES, FPS, DPI = P["frames"], P["fps"], P["dpi"]
FIGSIZE = (P["size"], P["size"])
MAX_POINTS, SPHERE_N = P["points"], P["n"]

ELEV = 25
TURNS = 1.0

# ---------- NOUVEAU : taille des boules vert/rouge ----------
# fraction du rayon du corps (0.1 = un dixième, 0.01 = un centième)
MARKER_FRAC = 0.01
MARKER_N = 12            # finesse des mini-sphères

# ---------- NOUVEAU : couleurs de la trajectoire ----------
COL_VISIBLE = "darkorange"   # partie visible (trait plein)
COL_HIDDEN = "dimgray"       # partie cachée derrière le corps (pointillés)

A = np.loadtxt('Geodesics.dat', ndmin=2)
B = np.loadtxt('dataforpy.txt')
r_s, r_moon, r_earth = B[0], B[1], B[2]

r, th, ph = A[:, 1], A[:, 2], A[:, 3]
x = r * np.sin(th) * np.cos(ph)
y = r * np.sin(th) * np.sin(ph)
z = r * np.cos(th)
r0 = r[0]

if r0 > r_earth:
    R_body, name = r_earth, "Earth"
else:
    R_body, name = r_s, "Schwarzschild radius"

show_rs = name != "Schwarzschild radius" and r_s > 1e-3 * r0
L = 1.1 * max(np.abs(x).max(), np.abs(y).max(), np.abs(z).max(), R_body)

step = max(1, len(x) // MAX_POINTS)
xs, ys, zs = x[::step], y[::step], z[::step]
idx = np.linspace(1, len(xs), N_FRAMES).astype(int)

R_marker = MARKER_FRAC * R_body


def sphere(R, n, c=(0, 0, 0)):
    u = np.linspace(0, 2*np.pi, n)
    v = np.linspace(0, np.pi, n)
    return (c[0] + R*np.outer(np.cos(u), np.sin(v)),
            c[1] + R*np.outer(np.sin(u), np.sin(v)),
            c[2] + R*np.outer(np.ones_like(u), np.cos(v)))


def draw_ball(ax, c, R, color, alpha=1.0, n=MARKER_N, zorder=5):
    X, Y, Z = sphere(R, n, c)
    return ax.plot_surface(X, Y, Z, color=color, alpha=alpha,
                           linewidth=0, shade=True, zorder=zorder)


def occluded(px, py, pz, elev, azim, R):
    """True pour les points masques par le corps vu depuis la camera."""
    e, a = np.radians(elev), np.radians(azim)
    d = np.array([np.cos(e)*np.cos(a), np.cos(e)*np.sin(a), np.sin(e)])
    depth = px*d[0] + py*d[1] + pz*d[2]
    perp2 = px**2 + py**2 + pz**2 - depth**2
    inside_disk = perp2 < R**2
    front_surface = np.sqrt(np.clip(R**2 - perp2, 0, None))
    return inside_disk & (depth < front_surface)


def split_visible(px, py, pz, elev, azim):
    m = occluded(px, py, pz, elev, azim, R_body)
    nan = np.nan
    vis = (np.where(m, nan, px), np.where(m, nan, py), np.where(m, nan, pz))
    hid = (np.where(m, px, nan), np.where(m, py, nan), np.where(m, pz, nan))
    return vis, hid


def make_scene(figsize, n):
    fig = plt.figure(figsize=figsize)
    # computed_zorder=False : on controle nous-memes l'ordre d'affichage
    ax = fig.add_subplot(projection='3d', computed_zorder=False)
    X, Y, Z = sphere(R_body, n)
    ax.plot_surface(X, Y, Z, color="lightsteelblue", alpha=0.55,
                    linewidth=0, antialiased=True, zorder=1)
    if show_rs:
        Xs, Ys, Zs = sphere(r_s, n)
        ax.plot_surface(Xs, Ys, Zs, color="black", alpha=0.8, linewidth=0, zorder=2)
    # boule verte (depart) : taille = MARKER_FRAC * rayon du corps
    draw_ball(ax, (x[0], y[0], z[0]), R_marker, "limegreen", zorder=6)
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_zlim(-L, L)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)"); ax.set_zlabel("z (m)")
    ax.set_title(f"Geodesic around {name}")
    return fig, ax


# ---------- rotating 3D animation ----------
fig, ax = make_scene(FIGSIZE, SPHERE_N)
line_vis, = ax.plot([], [], [], lw=2.5, color=COL_VISIBLE, zorder=4)
line_hid, = ax.plot([], [], [], lw=1.5, color=COL_HIDDEN, ls=(0, (2, 2)), zorder=3)
red_ball = [None]


def update(f):
    k = idx[f]
    azim = 30 + 360 * TURNS * f / N_FRAMES
    ax.view_init(elev=ELEV, azim=azim)

    px, py, pz = xs[:k], ys[:k], zs[:k]
    vis, hid = split_visible(px, py, pz, ELEV, azim)
    line_vis.set_data_3d(*vis)
    line_hid.set_data_3d(*hid)

    # boule rouge : meme taille relative, semi-transparente si derriere le corps
    if red_ball[0] is not None:
        red_ball[0].remove()
    c = (xs[k-1], ys[k-1], zs[k-1])
    hidden = occluded(np.array(c[0]), np.array(c[1]), np.array(c[2]),
                      ELEV, azim, R_body)
    red_ball[0] = draw_ball(ax, c, R_marker, "red",
                            alpha=0.35 if hidden else 1.0, zorder=7)
    return line_vis, line_hid


ani = animation.FuncAnimation(fig, update, frames=N_FRAMES, blit=False)

def progress(i, n):
    if i % 20 == 0:
        print(f"frame {i}/{n}", flush=True)

if animation.writers.is_available("ffmpeg"):
    writer = animation.FFMpegWriter(
        fps=FPS, codec="libx264",
        extra_args=["-crf", "18", "-preset", "medium", "-pix_fmt", "yuv420p"])
    ani.save("geodesic.mp4", writer=writer, dpi=DPI, progress_callback=progress)
    print("saved geodesic.mp4")
else:
    ani.save("geodesic.gif", writer=animation.PillowWriter(fps=FPS),
             dpi=DPI, progress_callback=progress)
    print("ffmpeg not found, saved geodesic.gif instead")
plt.close(fig)

# ---------- static plot ----------
fig, ax = make_scene((8, 8), SPHERE_N)
ax.view_init(elev=ELEV, azim=30)
vis, hid = split_visible(x, y, z, ELEV, 30)
ax.plot(*vis, lw=2, color=COL_VISIBLE, zorder=4, label="trajectory (visible)")
ax.plot(*hid, lw=1.2, color=COL_HIDDEN, ls=(0, (2, 2)), zorder=3,
        label="trajectory (derrière le corps)")
draw_ball(ax, (x[-1], y[-1], z[-1]), R_marker, "red", zorder=7)
ax.legend()
fig.savefig("static.png", dpi=150)
print("saved static.png")