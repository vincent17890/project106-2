from __future__ import division
import numpy as np
import networkx as nx
import math

def build_sbm(p, q, n, k):
  '''
  Build SBM adjacency matrix

  Parameters
  ----------
  p : float
      probability of edges when two nodes are in same community
  q : float
      probability of edges when two nodes are in difference communities
  n : integer
      number of nodes
  k : integer
      number of communities

  Return
  A : np.ndarray
      n x n adjacency matrix
  '''
  node_set = range(n)
  Z = np.zeros((1, n))
  A = np.zeros((n, n))
  
  # assign community evenly
  for node in node_set:
    for i in range(k):
      if node < (i+1)*int(n/k) and node >= i*int(n/k):
        Z[0, node] = i

  # build adjacency matrix
  for i in node_set:
    for j in range(i+1, n):
      if Z[0, i] == Z[0, j]:
        A[i, j] = np.random.binomial(1, p)
      else:
        A[i, j] = np.random.binomial(1, q)

  for i in node_set:
    for j in range(i):
      A[i, j] = A[j, i]

  return A

def output_graph(p, q, n, k):
  filename = 'sbm_' + str(int(p*10)) + '_' + str(int(math.log(n))) + '_' + str(k) + '.gml' 
  A = build_sbm(p, q, n, k)
  G = nx.from_numpy_matrix(A)
  nx.write_gml(G, filename)

def main():
 ps = [0.1, 0.5]
 qs = [0.01]
 ns = [100, 1000]
 ks = [2, 5]

 for p in ps:
  for q in qs:
    for n in ns:
      for k in ks:
        output_graph(p, q, n, k)
  
if __name__ == '__main__':
  main()
