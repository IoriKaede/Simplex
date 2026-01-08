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
def PrimalSimplex(A, b, c, B):
  n = len(A[0])
  N = list(range(n))
  for i in N.copy():
    if i in B:
      N.remove(i)
    else:
      pass
  count = 0
  while count < 9999:

    A_B = A[:, B]
    x = np.zeros(n)
    try:
      x[B] = np.linalg.solve(A_B, b)
      c_B = c[B]
      y = np.linalg.solve(A_B.T, c_B)
    except np.linalg.LinAlgError:
      return {"status": "singular matrix"}

    c_bar = np.zeros(n)
    c_bar[N] = c[N] - y.T @ A[:, N]

    if np.all(c_bar > -1e-7):
      return {"status": "optimal", "primal": x}

    k = min(j for j, v in enumerate(c_bar) if v < -1e-7)  # choosing the most negative k

    d = np.zeros(n)
    d[B] = np.linalg.solve(-A_B, A[:, k])
    d[k] = 1

    if np.all(d > -1e-7):
      return {"status": "unbounded", "primal": x}

    ratio = []  # ratio test
    for j in B:
      if d[j] < 0:
        ratio.append(-(x[j] / d[j]))
      else:
        ratio.append(float("inf"))
    l = np.argmin(ratio)

    N[N.index(k)] = B[l]
    count += 1
    B[l] = k

  return {"status" : "limit reached", "primal" : x }


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
  result = PrimalSimplex(A, b, c, list(basis))
  return result
  # So far we just print it.
  #print('Input LP:')
  #print(lp)


if __name__ == '__main__':

  file_name = 'orig-basis-blend.json'
  with open(file_name, 'r') as fp:
    lp = LP(fp.read())
    sol = solve(lp)


