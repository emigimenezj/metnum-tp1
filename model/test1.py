import numpy as np
from numpy.linalg import solve

a = np.array([[1, 0, -0.76, -6.84, 0],
              [0, 1, 0, 0, 0],
              [-2.28, 0, 1, -2.28, 0],
              [-2.28, 0, -1.52, 1, 0],
              [-2.28, 0, -0.76, -6.84, 1]])

b = np.ones(5)
print(a)
print(b)
sol = solve(a, b)

print(sol)
print("check")
print(a.dot(sol) - b)
normalized = sol / sum(sol)

print(normalized)
print(sum(normalized))
print(sum(np.array([-0.542893,
                    -2.59598,
                    0.588602,
                    0.452218,
                    3.09806])))
