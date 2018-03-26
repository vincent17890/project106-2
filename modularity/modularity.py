from __future__ import division
import numpy as np
import networkx as nx
from collections import deque
import 
'''
This program implement the modularity maximaization algorithm developed in
Modularity and community structures in networks by M. E. J. Newman
The binary version of this algorithm has been tested on karate dataset and reveal the true splitting
'''
def get_modularity_matrix(G):
  A = nx.to_numpy_matrix(G)
  n = A.shape[0]
  B = np.zeros((n, n))
  
  I = range(n)
  J = range(n)

  # the m is here already twice the number of edges
  m = sum([x[1] for x in nx.degree(G)])
  for i in I:
    for j in J:
      # Since the node number start from 1, we need i+1 and j+1
      ki = nx.degree(G, i+1)
      kj = nx.degree(G, j+1)
      B[i, j] = A[i, j] - ki*kj/m 
  return B

def get_submodularity_matrix(G, B=None, node_list=None):
  if node_list == None:
    return get_modularity_matrix(G)
  if B is None:
    B = get_modularity_matrix(G)
 
  # get indices of node in node_list
  indices = [list(G).index(node) for node in node_list]

  B_node_list = B[indices, :][:, indices]
  ng = len(indices)
  B_g = np.zeros((ng, ng), dtype=float)
  
  for i in range(ng):
    for j in range(ng):
      if i == j:
        B_g[i, j] = B_node_list[i, j] - sum(B_node_list[i])
      else:
        B_g[i, j] = B_node_list[i, j]

  return B_g

def get_delta_Q(B_g, a):
  '''
  Calculate delta modularity of a partition
  
  Parameters
  ----------
  B_g : np.matrix 
        submodularity matrix
  a : np.matrix
      community assignment
  
  Return
  ------
  dQ : float
       delta Q
  '''
  dQ = (a.dot(B_g)).dot(a.T)
  return dQ

def cmty_detect_binary(B):
  n = B.shape[0]
  cmty = np.zeros(n)
  
  w, v = np.linalg.eig(B)
  N = np.argmax(w)
  print w[N]
  print v[:, N]

  for i in range(n):
    if v[i, N] > 0:
      cmty[i] = 1
    else:
      cmty[i] = -1
  print cmty
  return cmty.tolist()

def divide(G, community_dict, community_index, B):
  '''
  Bisection of community in G

  Parameters
  ----------
  G : nx.Graph
      The network of interest
  
  Returns
  -------
  tuple
    If the community is indivisible, return (None, None)
    If the community is divisible, return 
    (node list of the first subgroup, node list of the original network)

  '''
  node_list = tuple(x for x in community_dict if community_dict[x] == community_index)
  B_g = get_submodularity_matrix(G, B, node_list)
  
  w, v = np.linalg.eig(B_g)
  N = np.argmax(w)
  top_eig_vec = v[:, N]
  top_eig_val = w[N]

  if top_eig_val > 0:
    # divisible
    a = np.asmatrix([1 if vi > 0 else -1 for vi in top_eig_vec])
    dQ = get_delta_Q(B_g, a)
    
    g1_node = np.array([node_list[i] for i in range(len(node_list)) if a[0, i] == 1])
    if len(g1_node) == 0 or len(g1_node) == len(node_list):
      return (None, None)
    
    if dQ > 0:
      return (g1_node, node_list)

    return (None, None)

  return (None, None)

def partition(G):
  ## Preprocessing
  G = nx.convert_node_labels_to_integers(G, first_label=1, label_attribute='node_name')
  node_name = nx.get_node_attributes(G, 'node_name')

  ## Only support unweighted network
  nx.set_edge_attributes(G=G, name='weight', values={edge:1 for edge in G.edges})

  B = get_modularity_matrix(G)
 
  ## Set flags for divisibilities of communities
  ## Inital community is divisibl
  divisible_community = deque([0])

  ## All nodes are in the same community
  community_dict = {node:0 for node in G}

  community_counter = 0
  while(len(divisible_community) > 0):
    community_index = divisible_community.popleft()
    g1, full_community = divide(G, community_dict, community_index, B)
    if g1 is None:
      # indisible, continue to next
      continue
    
    g2 = set(full_community).difference(set(g1))
    
    community_counter += 1
    divisible_community.append(community_counter)
    for u in g1:
      community_dict[u] = community_counter
    
    community_counter += 1
    divisible_community.append(community_counter)
    for u in g2:
      community_dict[u] = community_counter
  
  return {node_name[u]: community_dict[u] for u in G} 

def main():
  graphpath = './polblogs/polblogs.gml'
  #G = nx.read_edgelist(graphpath)
  G = nx.read_gml(graphpath, label='id')
  
  cmty = partition(G)
  community_set = set([cmty[x] for x in cmty])
  print len(community_set)
if __name__ == '__main__':
  main()
