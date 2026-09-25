import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

A = np.loadtxt('Geodesics.txt')
B = np.loadtxt('dataforpy.dat')   # not used yet — plot it too if it's meant for comparison

x, y, z = (A[:,1]*np.sin(A[:,2])*np.cos(A[:,3]),
           A[:,1]*np.sin(A[:,2])*np.sin(A[:,3]),
           A[:,1]*np.cos(A[:,2]))

fig, ax = plt.subplots()
ax.plot(x, y, lw=1, label="trajectory")

r_s = 8.87e-3  # or compute it from G, M, c as before
ax.add_patch(Circle((0, 0), r_s, fill=False, color="black", label="horizon"))

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.axis("equal")
ax.set_title("Geodesic trajectory")
ax.legend()
plt.show()