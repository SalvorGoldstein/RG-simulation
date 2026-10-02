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

def sphere(R, n):
    u = np.linspace(0, 2*np.pi, n)
    v = np.linspace(0, np.pi, n)
    return (R*np.outer(np.cos(u), np.sin(v)),
            R*np.outer(np.sin(u), np.sin(v)),
            R*np.outer(np.ones_like(u), np.cos(v)))

def make_scene(figsize, n):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(projection='3d')
    X, Y, Z = sphere(R_body, n)
    ax.plot_surface(X, Y, Z, color="steelblue", alpha=0.6, linewidth=0, antialiased=True)
    if show_rs:
        Xs, Ys, Zs = sphere(r_s, n)
        ax.plot_surface(Xs, Ys, Zs, color="black", alpha=0.8, linewidth=0)
    ax.scatter([x[0]], [y[0]], [z[0]], color="green", s=35)
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_zlim(-L, L)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)"); ax.set_zlabel("z (m)")
    ax.set_title(f"Geodesic around {name}")
    return fig, ax

# ---------- rotating 3D animation ----------
fig, ax = make_scene(FIGSIZE, SPHERE_N)
line, = ax.plot([], [], [], lw=1.5, color="C0")
point, = ax.plot([], [], [], "o", color="red", ms=6)

def update(f):
    k = idx[f]
    line.set_data_3d(xs[:k], ys[:k], zs[:k])
    point.set_data_3d([xs[k-1]], [ys[k-1]], [zs[k-1]])
    ax.view_init(elev=ELEV, azim=30 + 360 * TURNS * f / N_FRAMES)
    return line, point

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
ax.plot(x, y, z, lw=1, label="trajectory")
ax.legend()
fig.savefig("static.png", dpi=150)
print("saved static.png")