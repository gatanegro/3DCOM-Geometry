import numpy as np
from scipy.spatial.distance import pdist, squareform

# --- Shell 3 Data ---
shell3_primes = [31, 37, 41, 43, 47, 53, 59, 61,
                 67, 71, 73, 79, 83, 89, 97, 101, 103, 107]
shell3_roots = [4, 1, 5, 7, 2, 8, 5, 7, 4, 8, 1, 7, 2, 8, 7, 2, 4, 8]
LZ = 1.23498228


def map_prime_to_sphere(p, root, LZ_attr):
    theta = (root / 9) * 2 * np.pi
    remainder = p % LZ_attr
    normalized_rem = remainder / LZ_attr
    phi = np.arcsin(2 * normalized_rem - 1)
    return theta, phi


# Generate coordinates for Shell 3 points
R = 3  # Shell radius
xs, ys, zs = [], [], []
for p, root in zip(shell3_primes, shell3_roots):
    theta, phi = map_prime_to_sphere(p, root, LZ)
    x = R * np.cos(phi) * np.cos(theta)
    y = R * np.cos(phi) * np.sin(theta)
    z = R * np.sin(phi)
    xs.append(x)
    ys.append(y)
    zs.append(z)

points = np.vstack((xs, ys, zs)).T

# Normalize points to direction vectors
norm_points = points / np.linalg.norm(points, axis=1, keepdims=True)

# Compute angles between all pairs
cosine_matrix = np.dot(norm_points, norm_points.T)
angle_matrix = np.arccos(np.clip(cosine_matrix, -1, 1)) * 180 / np.pi

# Characteristic dodecahedron angles
angle_adjacent = 116.565  # Angle between adjacent vertices
angle_face = 41.81        # Angle across a pentagonal face

# Check for angles close to these values
print("Checking for dodecahedron angles in Shell 3 (≈116.565° and ≈41.81°):")
adjacent_count = 0
face_count = 0

for i in range(len(angle_matrix)):
    for j in range(i+1, len(angle_matrix)):
        angle = angle_matrix[i, j]
        if 114 < angle < 118:  # Check for adjacent vertex angle
            adjacent_count += 1
            print(
                f"Adjacent angle found: {angle:.3f}° between {shell3_primes[i]} and {shell3_primes[j]}")
        if 40 < angle < 43:    # Check for face angle
            face_count += 1
            print(
                f"Face angle found: {angle:.3f}° between {shell3_primes[i]} and {shell3_primes[j]}")

print(f"\nTotal adjacent (~116.565°) angles found: {adjacent_count}")
print(f"Total face (~41.81°) angles found: {face_count}")

# In a dodecahedron, each of the 20 vertices has 3 adjacent vertices, so there are (20*3)/2 = 30 edges.
# We expect to find many angles near 116.565°.
