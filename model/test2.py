5
12
1
5
4
1
3
4
1
4
2
1
4
3
3
1
2
4
1
3
2
5
4
5
3
5

import numpy as np

w = np.zeros((5, 5))

for i,j in [(1, 5),
          (4, 1),
          (3, 4),
          (1, 4),
          (2, 1),
          (4, 3),
          (3, 1),
          (2, 4),
          (1, 3),
          (2, 5),
          (4, 5),
          (3, 5)]:

    w[j-1][i-1]=1.0

EPS = 1e-8
d = np.diag(np.ones(5))
for i in range(0, 5):
    c = np.sum(w[:,i])
    d[i][i] = 0. if np.abs(c) < EPS else 1.0/c

p = 0.76
I = np.diag(np.ones(5))
print("P")
print(p)
print("W")
print(w)
print("D")
print(d)
print("I")
print(I)
print("Cuenta")
print(p * np.matmul(w, d))

a = I - (p * np.matmul(w, d))

print("Cuenta final")
print(a)


from numpy.linalg import solve

b = np.ones(5)

sol = solve(a, b)
print("SOL")
print(sol)
print("SOL suma")
print(np.sum(sol))
normalized = sol/np.sum(sol)
print("CHEK")
print(np.sum(normalized))

print("Solucion normalizada")

print(normalized)
