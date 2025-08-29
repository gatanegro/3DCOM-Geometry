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
    return np.log(prime) / np.log(LZ_attr)


# --- Generate Primes and Assign to Energy Bands ---
primes = list(primerange(2, 100))
energy_bands = {}
band_width = 2.0  # Group layers into bands of this width

for p in primes:
    L = assign_layer_log(p, LZ)
    band_id = int(L // band_width)  # Group by integer division
    if band_id not in energy_bands:
        energy_bands[band_id] = []
    energy_bands[band_id].append(p)

# Print the new energy band assignment
print("Primes assigned to Energy Bands (Band Width = 2.0):")
for band, primes_in_band in sorted(energy_bands.items()):
    print(f"Band {band}: {primes_in_band}")

# --- Map Primes on Their Energy Band Sphere ---


def map_prime_to_sphere(p, LZ_attr):
    root = digital_root(p)
    theta = (root / 9) * 2 * np.pi  # Longitude

    remainder = p % LZ_attr
    normalized_rem = remainder / LZ_attr
    phi = np.arcsin(2 * normalized_rem - 1)  # Latitude

    return theta, phi


# --- 3D Plot for Each Energy Band ---
for band_id, primes_in_band in energy_bands.items():
    if len(primes_in_band) < 3:
        continue

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    R = band_id  # Use band ID as the radius for visualization

    xs, ys, zs = [], [], []
    for p in primes_in_band:
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

    # Find and plot the convex hull
    points = np.vstack((xs, ys, zs)).T
    try:
        hull = ConvexHull(points)
        for simplex in hull.simplices:
            simplex = np.append(simplex, simplex[0])
            ax.plot(points[simplex, 0], points[simplex, 1],
                    points[simplex, 2], 'r-', linewidth=2, alpha=0.5)
        plt.title(
            f'Energy Band {band_id}\nPrimes: {primes_in_band}\nShape: {len(hull.vertices)}-sided Polygon')
    except:
        plt.title(
            f'Energy Band {band_id}\nPrimes: {primes_in_band}\n(Points are co-linear)')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    max_range = R * 1.5
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)
    ax.set_zlim(-max_range, max_range)
    plt.show()
