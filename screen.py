import numpy as np
import matplotlib.pyplot as plt

A = np.loadtxt('Geodesics.dat', ndmin=2)
B = np.loadtxt('dataforpy.txt')
r_s, r_moon, r_earth = B[0], B[1], B[2]

r, th, ph = A[:, 1], A[:, 2], A[:, 3]
x = r * np.sin(th) * np.cos(ph)
y = r * np.sin(th) * np.sin(ph)
z = r * np.cos(th)

r0 = r[0]

# Pick the central body that matches the scale of r0:
# the largest known radius that is still smaller than r0.
if r0 > r_earth:
    R_body, name = r_earth, "Earth"
else:
    R_body, name = r_s, "Schwarzschild radius"

def sphere(R, n=40):
    u = np.linspace(0, 2*np.pi, n)
    v = np.linspace(0, np.pi, n)
    return (R*np.outer(np.cos(u), np.sin(v)),
            R*np.outer(np.sin(u), np.sin(v)),
            R*np.outer(np.ones_like(u), np.cos(v)))

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(projection='3d')

ax.plot(x, y, z, lw=1, label="trajectory")
ax.scatter(x[0], y[0], z[0], color="green", s=30, label=f"start (r0 = {r0/R_body:.2f} R)")

X, Y, Z = sphere(R_body)
ax.plot_surface(X, Y, Z, color="steelblue", alpha=0.6, linewidth=0)

# Also show r_s if it is a visible fraction of the picture
if name != "Schwarzschild radius" and r_s > 1e-3 * r0:
    Xs, Ys, Zs = sphere(r_s)
    ax.plot_surface(Xs, Ys, Zs, color="black", alpha=0.8, linewidth=0)

# Equal axes, scaled to the data and the body
L = 1.1 * max(np.abs(x).max(), np.abs(y).max(), np.abs(z).max(), R_body)
ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_zlim(-L, L)
ax.set_box_aspect((1, 1, 1))

ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_zlabel("z (m)")
ax.set_title(f"Geodesic trajectory around {name}")
ax.legend()
plt.show()