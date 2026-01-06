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
  for con in lp.constraints:
    b = np.array(con['rhs'])
    c = np.array(lp.objective)
  print(m,n,A,b,c)

  #while limit:
    n = list(range((lp.num_columns)+1))
    n = n.remove(0)
    basis = lp.basis
    print(n)
    for i in n:
      if i in basis:
        n.remove(i)
      else:
        pass

  # So far we just print it.
  print('Input LP:')
  print(lp)


if __name__ == '__main__':

  file_name = 'BT-Example-3.5-std.json'
  with open(file_name, 'r') as fp:
    lp = LP(fp.read())
    solve(lp)

