import numpy as np
from sympy import primerange

# --- Function to compute digital root ---


def digital_root(n):
    return (n - 1) % 9 + 1


# Generate first 100 primes
primes = list(primerange(2, 542))[:100]  # First 100 primes

# Atomic shell capacities (electrons per shell)
shell_capacities = [2, 8, 18, 32, 50]  # K, L, M, N, O shells

# Assign primes to atomic-like shells
shells = {}
current_index = 0
for shell_num, capacity in enumerate(shell_capacities, 1):
    shell_primes = primes[current_index:current_index + capacity]
    shells[shell_num] = shell_primes
    current_index += capacity
    if current_index >= len(primes):
        break

# Print the results
print("Atomic-like Shells for Prime Numbers:")
for shell_num, primes_in_shell in shells.items():
    print(
        f"Shell {shell_num} (Capacity: {shell_capacities[shell_num-1]}): {primes_in_shell}")
    print(f"  Number of primes: {len(primes_in_shell)}")
    digital_roots = [digital_root(p) for p in primes_in_shell]
    print(f"  Digital roots: {digital_roots}")
    print()
