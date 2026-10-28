import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

base = "/home/adam/Documents/study/prog/moredim/RG-simulation/python/"
data = np.loadtxt(base + "dataforpy.dat")
n = int(data[4])

for i in range(n):
    Geodesic = np.loadtxt(base + "data/Geodesics" + str(i) + ".dat")
    r, theta, phi = Geodesic[:, 1], Geodesic[:, 2], Geodesic[:, 3]
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    plt.plot(x, y)

axes = plt.gca()
axes.add_patch(patches.Circle((0, 0), radius=data[0], color='k', zorder=2))
axes.set_aspect('equal', adjustable='box')

# new: extra dimension w (column 4) against t (column 0)
plt.figure()
for i in range(n):
    Geodesic = np.loadtxt(base + "data/Geodesics" + str(i) + ".dat")
    plt.plot(Geodesic[:, 0], Geodesic[:, 4])
plt.xlabel("x0 (t)")
plt.ylabel("x4 (w)")

plt.show()