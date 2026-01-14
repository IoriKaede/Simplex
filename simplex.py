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

from lp import LP
import numpy as np


def PrimalSimplex(A, b, c, B):
  # Set N:= {1,2, ...,k} \ B
  n = len(A[0])
  N = list(range(n))
  for i in N.copy():
    if i in B:
      N.remove(i)
    else:
      pass
  count = 0
  while count < 9999:
    print(count)
    #Solve A_{*,B} x_B = b for x_B ∈ R^B and set x_N := O_N
    #Solve y⊺A_{*,B} = c_B for y ∈ R^m and compute ¯c_j := c_j − y⊺A_{*,B} for all j ∈ N
    A_B = A[:, B]
    x = np.zeros(n)
    try:
      x[B] = np.linalg.solve(A_B, b)
      c_B = c[B]
      y = np.linalg.solve(A_B.T, c_B) # y is the optimal dual solution
    except np.linalg.LinAlgError:
      return {"status": "singular matrix"}

    c_bar = np.zeros(n)
    c_bar[N] = c[N] - y.T @ A[:, N]
    #If ¯c_N ≥ O then return (“optimal”, solution x, optimal basis B)
    #if np.all(c_bar > -1e-7):
    if np.all(c_bar[N] > -1e-7):
        return {"status": "optimal", "primal": x, "dual" : y, "basis":B}  #add basis as return
                    # small change: c_bar[N] instead of c_bar
                    # In c_bar we also have c_bar[B] which should be zero's, but due to noise or if indexing is mixed up, it can get <0

    #Choose the most negative k from {j ∈ N | ¯c_j < 0}
    # k = min(j for j, v in enumerate(c_bar) if v < -1e-7)
    k = min(j for j in N if c_bar[j] < -1e-7)
                    # small change: We shouldn't risk having a j from B.

    #Solve A_{*,B} d_B = A_{*,j} for d_B ∈ R^B and set d_k := 1, leave d_{N\{k}} := 0.
    d = np.zeros(n)
    d[B] = np.linalg.solve(-A_B, A[:, k])
    d[k] = 1

    #If d ≥ O_N then return (“unbounded“, solution x, unbounded direction d).
    if np.all(d > -1e-7):
      return {"status": "unbounded", "primal": x, "dual" : y, "ray":d, "basis":B} #direction included

    #Compute θ* := min{−x_j/d_j | j ∈ B, d_j < 0}
    ratio = []  # ratio test
    for j in B:
      # if d[j] < 0:
      if d[j] < -1e-7:
        ratio.append(-(x[j] / d[j]))
      else:
        ratio.append(float("inf"))
    #Choose l from {j ∈ B | xj = −dj · θ*}
    l = np.argmin(ratio)

    #Update B := (B \ {ℓ}) ∪ {k}
    # N[N.index(k)] = B[l]
    # B[l] = k
    B_l = B[l]
    B[l] = k
    N.remove(k)
    N.append(B_l)
                    # small change: Changed N[N.index(k)], since index doesn't matter. Have to make sure it actually removes k.
    count += 1
  return {"status" : "limit reached", "primal" : x , "dual" : y, "basis":B}

def PrimalSimplex_without_basis(A,b,c):
  m = lp.num_rows
  n = lp.num_columns
  A_bar = np.hstack((A , np.identity(m)))
  b_bar = b
  c_bar = np.zeros(n+m)
  c_bar[n:] = 1
  B_bar = []
  for i in range(n,n+m):
    B_bar.append(i)
  print("phase 1")
  dic = PrimalSimplex(A_bar, b_bar, c_bar, B_bar)
  print("after phase 1")
  status = dic["status"]
  x_prime = dic["primal"]
  B_prime = dic["basis"]
  B_p = list(B_prime)
  print(x_prime)

  if status =="limit reached":
    return {"status":"limit reached"}
  if c_bar.T @ x_prime > 0 :
    return {"status":"infeasible"}






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
  for i in range(len(b)):
    if b[i] <0:
      A[i] = - A[i]
      b[i] = - b[i]
            # we use arrays so we can use matrix multiplication
  # basis = lp.basis
            #small change: might not work if there is no basis given
  try:
    basis = lp.basis
  except KeyError:
    basis = None
  if basis is None:
    basis = PrimalSimplex_without_basis(A, b, c)
    if basis is None:
      return {"status": "infeasible"}

            #small change: Use slack variables if no basis provided


  result = PrimalSimplex(A, b, c, list(basis))
  return result
  # # So far we just print it.
  # print('Input LP:')
  # print(lp)


if __name__ == '__main__':

  file_name = 'orig-std-adlittle.json'
  with open(file_name, 'r') as fp:
    lp = LP(fp.read())
    sol = solve(lp)
    # print(lp.primal_value(sol["primal"]))
    print(sol)