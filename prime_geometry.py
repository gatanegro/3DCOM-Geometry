import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull

# --- 3DCOM Constants ---
LZ = 1.23498228

# --- Functions ---


def digital_root(n):
    return (n - 1) % 9 + 1


def map_prime_to_sphere(p, LZ_attr):
    root = digital_root(p)
    theta = (root / 9) * 2 * np.pi  # Longitude
    remainder = p % LZ_attr
    normalized_rem = remainder / LZ_attr
    phi = np.arcsin(2 * normalized_rem - 1)  # Latitude
    return theta, phi


# --- Define the energy bands and their primes ---
bands = {
    6: [13, 17, 19],
    8: [31, 37, 41, 43],
    9: [47, 53, 59, 61, 67],
    10: [71, 73, 79, 83, 89, 97]
}

# --- 3D Plot for Each Energy Band ---
for band_id, primes_in_band in bands.items():
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    R = band_id  # Use band ID as the radius

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
                fontsize=12, ha='center', va='center')

    ax.scatter(xs, ys, zs, s=100, c='blue', alpha=0.8)

    # Find and plot the convex hull
    points = np.vstack((xs, ys, zs)).T
    try:
        hull = ConvexHull(points)
        for simplex in hull.simplices:
            simplex = np.append(simplex, simplex[0])
            ax.plot(points[simplex, 0], points[simplex, 1],
                    points[simplex, 2], 'r-', linewidth=3, alpha=0.6)
        plt.title(
            f'Energy Band {band_id}\nPrimes: {primes_in_band}\nShape: {len(hull.vertices)}-sided Polygon', fontsize=14)
    except Exception as e:
        plt.title(
            f'Energy Band {band_id}\nPrimes: {primes_in_band}\n(No convex hull: {e})')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    max_range = R * 1.5
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)
    ax.set_zlim(-max_range, max_range)
    plt.show()
