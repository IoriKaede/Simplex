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
  basis = lp.basis
  B = list(basis)
  N = list(range(n))

  count = 0
  while count < 9999:
    for i in N.copy():
      if i in basis:
        N.remove(i)
      else:
        pass

    A_B = A[:, B]
    x_B = np.linalg.solve(A_B, b)
    c_B = c[B]
    y = np.linalg.solve(A_B.T, c_B.T)

    x_N = np.zeros(len(N))
    x = np.append(x_B, x_N)

    c_j = np.array([])
    for j in N:
      c_j= np.append(c_j , c[j] - y.T @ A[:, j])

    if np.all(c_j >= 0):
      return "optimal", x, basis

    k = np.argmin(c_j)  #choosing the most negative k

    d_B = np.linalg.solve(-A_B, A[:, k])
    d = np.zeros(n)
    d[B] = d_B
    d[k] = 1

    if np.all(d >= 0):
      return "unbounded", x, d

    ratio = np.array([])      #ratio test
    for j in B:
      if d[j] < 0:
        ratio = np.append(ratio, -(x[j] / d[j]))
      else:
        pass
    #theta_star = np.min(ratio)

    l = np.argmin(ratio)
    B[l] = k
    count += 1
  return "limit reached", x, B
  # So far we just print it.
  #print('Input LP:')
  #print(lp)


if __name__ == '__main__':

  file_name = 'BT-Example-3.5-std.json'
  with open(file_name, 'r') as fp:
    lp = LP(fp.read())
    solve(lp)

