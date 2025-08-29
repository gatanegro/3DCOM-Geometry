import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sympy import primerange
from scipy.spatial import ConvexHull

# --- 3DCOM Constants ---
LZ = 1.23498228

# --- Functions ---


def digital_root(n):
    return (n - 1) % 9 + 1


def assign_layer_log(prime, LZ_attr):
    """Assign a prime to a layer based on the logarithm of its value."""
    # Layer = log_base(LZ) of (prime) = ln(prime) / ln(LZ)
    return np.log(prime) / np.log(LZ_attr)


# --- Generate Primes and Assign Layers (Logarithmic) ---
primes = list(primerange(2, 100))  # Primes up to 100
layers_log = {}
for p in primes:
    L = assign_layer_log(p, LZ)
    # We'll group them into integer bins for visualization
    L_bin = int(np.floor(L))
    if L_bin not in layers_log:
        layers_log[L_bin] = []
    layers_log[L_bin].append(p)

# Print the new layer assignment
print("Primes assigned to layers (based on ln(p) / ln(LZ)):")
for L_bin, primes_in_L in sorted(layers_log.items()):
    print(f"Layer ~{L_bin}: {primes_in_L}")

# --- Map Primes on their Layer Sphere ---


def map_prime_to_sphere(p, LZ_attr):
    root = digital_root(p)
    theta = (root / 9) * 2 * np.pi  # Longitude (0 to 2π)

    remainder = p % LZ_attr
    normalized_rem = remainder / LZ_attr  # Map to [0, 1]
    # Map [0,1] to [-π/2, π/2] for Latitude
    phi = np.arcsin(2 * normalized_rem - 1)

    return theta, phi


# --- 3D Plot for Dense Layers ---
for L_bin, primes_in_L in layers_log.items():
    if len(primes_in_L) < 3:  # Only plot layers with enough primes to form a shape
        continue

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    R = L_bin  # Radius of the sphere for this layer bin

    xs, ys, zs = [], [], []
    for p in primes_in_L:
        theta, phi = map_prime_to_sphere(p, LZ)
        x = R * np.cos(phi) * np.cos(theta)
        y = R * np.cos(phi) * np.sin(theta)
        z = R * np.sin(phi)

        xs.append(x)
        ys.append(y)
        zs.append(z)
        ax.text(x, y, z, str(p), color='darkred',
                fontsize=10, ha='center', va='center')

    ax.scatter(xs, ys, zs, s=80, c='blue', alpha=0.8)

    # Find and plot the convex hull to see the geometric shape
    points = np.vstack((xs, ys, zs)).T
    try:
        hull = ConvexHull(points)
        for simplex in hull.simplices:
            simplex = np.append(simplex, simplex[0])  # close the polygon
            ax.plot(points[simplex, 0], points[simplex, 1],
                    points[simplex, 2], 'r-', linewidth=2, alpha=0.5)
        plt.title(
            f'Layer ~{L_bin}\nPrimes: {primes_in_L}\nShape: {hull.vertices.size}-sided Polygon')
    except:
        plt.title(
            f'Layer ~{L_bin}\nPrimes: {primes_in_L}\n(Points are co-linear)')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    max_range = R * 1.2
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)
    ax.set_zlim(-max_range, max_range)
    plt.show()
