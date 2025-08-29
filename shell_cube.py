import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull

# --- Shell 2 Data ---
shell2_primes = [5, 7, 11, 13, 17, 19, 23, 29]
shell2_roots = [5, 7, 2, 4, 8, 1, 5, 2]
LZ = 1.23498228


def map_prime_to_sphere(p, root, LZ_attr):
    # Use the digital root for longitude
    theta = (root / 9) * 2 * np.pi  # Longitude based on digital root
    # Use the prime value modulo LZ for latitude
    remainder = p % LZ_attr
    normalized_rem = remainder / LZ_attr
    phi = np.arcsin(2 * normalized_rem - 1)  # Map to latitude [-π/2, π/2]
    return theta, phi


# Plot Shell 2
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
R = 2  # Shell radius

xs, ys, zs = [], [], []
for p, root in zip(shell2_primes, shell2_roots):
    theta, phi = map_prime_to_sphere(p, root, LZ)
    x = R * np.cos(phi) * np.cos(theta)
    y = R * np.cos(phi) * np.sin(theta)
    z = R * np.sin(phi)
    xs.append(x)
    ys.append(y)
    zs.append(z)
    ax.text(x, y, z, f'{p}({root})', color='darkred',
            fontsize=10, ha='center', va='center')

ax.scatter(xs, ys, zs, s=100, c='blue', alpha=0.8)

# Find and plot the convex hull (should be a cube)
points = np.vstack((xs, ys, zs)).T
hull = ConvexHull(points)
for simplex in hull.simplices:
    simplex = np.append(simplex, simplex[0])
    ax.plot(points[simplex, 0], points[simplex, 1],
            points[simplex, 2], 'r-', linewidth=2, alpha=0.6)

plt.title('Shell 2 (L-Shell) Geometry: 8 Primes Forming a Cube\nPrime(Digital Root)')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
max_range = R * 1.5
ax.set_xlim(-max_range, max_range)
ax.set_ylim(-max_range, max_range)
ax.set_zlim(-max_range, max_range)
plt.show()

# Check if the shape is a cube by analyzing the hull
print(f"Shell 2 Convex Hull has {len(hull.vertices)} vertices.")
# A cube has 8 vertices and 12 edges
print(f"Number of edges (simplices): {len(hull.simplices)}")
