# ‐-‐-‐‐‐-‐‐-‐‐---‐-‐-‐---‐----‐----‐‐---‐‐---‐‐‐‐-‐----‐‐‐--‐‐‐--‐‐‐‐-‐‐‐‐-‐‐‐-‐-‐-‐‐‐--‐--‐---‐-‐-
#
# Implement the solve-function as specified in the assignment description. Feel free to use Linear
# Algebra routines from the numpy and scipy packages. However, usage of routines that actually solve
# LPs is not allowed.
#
# ‐‐-‐---‐‐‐-‐‐‐‐‐--‐‐-‐‐‐‐--‐-‐--‐-------‐-‐‐-‐‐-‐-‐--‐--‐-‐‐‐‐‐‐‐-‐‐-‐-‐‐----‐--‐‐-----‐---‐‐‐-‐--

'''{
  "sense": "minimize",
  "objective": [ -10.0, -12.0, -12.0, 0.0, 0.0, 0.0 ],
  "signs": [ 1, 1, 1, 1, 1, 1 ],
  "constraints": [
    { "coefficients": { "0": 1.0, "1": 2.0, "2": 2.0, "3": 1.0 }, "relation": "=", "rhs": 20 },
    { "coefficients": { "0": 2.0, "1": 1.0, "2": 2.0, "4": 1.0 }, "relation": "=", "rhs": 20 },
    { "coefficients": { "0": 2.0, "1": 2.0, "2": 1.0, "5": 1.0 }, "relation": "=", "rhs": 20 }
  ],
  "basis": [ 3, 4, 5 ]
}'''


import numpy as np
from lp import LP

def solve(lp):
  m = lp.num_rows
  n = lp.num_columns
  A = np.zeros((m, n))
  for i, con in enumerate(lp.constraints):
    for j, val in con['coefficients'].items():
      A[i, j] = val
  b_list = []
  for con in lp.constraints:
    b_list.append(float(con['rhs']))
  b = np.array(b_list)
  c = np.array(lp.objective)

  #while limit:
  N = list(range(n))
  print(N)
  basis = lp.basis
  for i in N.copy():
    if i in basis:
      N.remove(i)
    else:
      pass
  print(N)

  A_B = A[:, list(basis)]
  x_B = np.linalg.solve(A_B, b)
  c_B = c[list(basis)]
  y = np.linalg.solve(A_B.T, c_B.T)
  print(y)

  x_N = np.zeros(len(N))
  x = np.append(x_B, x_N)
  print(x)

  c_j = np.array([])
  for j in N:
    c_j= np.append(c_j , c[j] - y.T @ A[:, j])
  print(c_j)

  # ↓ not finished
  if c_j >= np.zeros(len(c_j)):
    print(x,B)
    return "optimal", x, basis

  # So far we just print it.
  print('Input LP:')
  print(lp)


if __name__ == '__main__':

  file_name = 'BT-Example-3.5-std.json'
  with open(file_name, 'r') as fp:
    lp = LP(fp.read())
    solve(lp)

