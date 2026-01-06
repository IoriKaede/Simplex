# ‐-‐-‐‐‐-‐‐-‐‐---‐-‐-‐---‐----‐----‐‐---‐‐---‐‐‐‐-‐----‐‐‐--‐‐‐--‐‐‐‐-‐‐‐‐-‐‐‐-‐-‐-‐‐‐--‐--‐---‐-‐-
#
# Implement the solve-function as specified in the assignment description. Feel free to use Linear
# Algebra routines from the numpy and scipy packages. However, usage of routines that actually solve
# LPs is not allowed.
#
# ‐‐-‐---‐‐‐-‐‐‐‐‐--‐‐-‐‐‐‐--‐-‐--‐-------‐-‐‐-‐‐-‐-‐--‐--‐-‐‐‐‐‐‐‐-‐‐-‐-‐‐----‐--‐‐-----‐---‐‐‐-‐--

import sys
import numpy
from lp import LP

def solve(lp):
  # So far we just print it.
  print('Input LP:')
  print(lp)


if __name__ == '__main__':

  if len(sys.argv) >= 2:
    # Read the file name from the command line (in case you call it like that:
    # python3 simplex.py orig-basis-blend.json
    file_name = sys.argv[1]
  else:
    # Or manually set a name here.
    file_name = 'orig-basis-blend.json'

  # First argument is the file name.
  file_name = sys.argv[1]
  with open(file_name, 'r') as fp:
    lp = LP(fp.read())
    solve(lp)

