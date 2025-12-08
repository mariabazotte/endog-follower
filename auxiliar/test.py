import numpy as np
from PolyRound.mutable_classes.polytope import Polytope
from PolyRound.api import PolyRoundApi
from PolyRound.settings import PolyRoundSettings

A = np.array([
    [ 1,  0],
    [ 0,  1],
    [-1,  0],
    [ 0, -1]
], dtype=np.float64)

b = np.array([1, 1, 0, 0], dtype=np.float64)

# -------------------------------
# Create Polytope
# -------------------------------
poly = Polytope(A, b)

# -------------------------------
# Transform and round polytope
# -------------------------------
# poly = PolyRoundApi.simplify_polytope(poly)
# poly = PolyRoundApi.transform_polytope(poly)
poly_rounded = PolyRoundApi.round_polytope(poly)

# -------------------------------
# Extract transformed/rounded results
# -------------------------------
A_new = np.array(poly_rounded.A.values)
b_new = np.array(poly_rounded.b.values)
S     = np.array(poly_rounded.transformation.values)
t     = np.array(poly_rounded.shift)

print("\nTransformed & Rounded A:\n", A_new)
print("Transformed & Rounded b:\n", b_new)
print("Transformation matrix S:\n", S)
print("Shift vector t:\n", t)

# -------------------------------
# Optional: sanity check
# Check if A * S + t ≈ A_new (approx)
# -------------------------------
check = np.allclose(A @ S , A_new, atol=1e-12)
print("\nSanity check (A @ S ≈ A_new):", check)

# pick x feasible in original polytope
x = np.array([0.5, 0.5])

# transformed point
y = S @ x + t

print("Original constraints satisfied:", np.all(A @ x <= b))
print("Transformed constraints satisfied:", np.all(A_new @ y <= b_new))