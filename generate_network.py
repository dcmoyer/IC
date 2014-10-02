#!/usr/bin/python
#
#  IC
#
#

import matplotlib
matplotlib.use('Agg')

import networkx as nx
import numpy as np
import pickle
import matplotlib.pyplot as pl

def diff(a, b):
    b = set(b)
    return [aa for aa in a if aa not in b]

#generate_forest_fire_network
#   
#   n       The number of nodes to be added to the network.
#   p       The forward burning parameter
#   r       The backward burning parameter
#   file_name The output filename.
'''@inproceedings{leskovec2005graphs,
  title={Graphs over time: densification laws, shrinking diameters and possible explanations},
  author={Leskovec, Jure and Kleinberg, Jon and Faloutsos, Christos},
  booktitle={Proceedings of the eleventh ACM SIGKDD international conference on Knowledge discovery in data mining},
  pages={177--187},
  year={2005},
  organization={ACM}
}'''
def generate_forest_fire_network(n,p = 0.3,r = 0.3 ,file_name='test_graph.p',directory_path='./cascade_data/',reverse_flag=False):
    
    #create a graph with one node
    G = nx.DiGraph()
    G.add_node(1)
    
    #iteratively add the rest of the nodes
    for new_node_number in range(1,n):
        G.add_node(new_node_number)
        
        #select ambassador
        w = np.random.random_integers(0,new_node_number - 1, 1)
        w = int(w)
        G.add_edge(new_node_number,w)
        
        
        #initial list of targets
        list_of_targets = [w]
        
        #for each target
        for x in list_of_targets:
            #find the neighbors of the target that are not already neighbors of the new node
            target_neighbors = diff(G.neighbors(x),G.neighbors(new_node_number))
            if(len(target_neighbors) > 0):
                
                
                #randomly select a random number of them
                #uniform selection of a binomial number of them
                forward_burn_number = np.random.binomial(len(target_neighbors),p,1)
                forward_burn_set = np.random.choice(target_neighbors,forward_burn_number, replace=False)
                
                #add more neighbors, with possible return edges
                for i in forward_burn_set:
                    #outward edges
                    G.add_edge(new_node_number,i)
                    #now search the ones we just added
                    list_of_targets.append(i)
                #end
            #end
            
            L = G.neighbors(new_node_number)
            L.append(new_node_number)
            target_neighbors = diff(G.predecessors(x),L)
            if(len(target_neighbors) == 0):
                continue
            forward_burn_number = np.random.binomial(len(target_neighbors),p*r,1)
            forward_burn_set = np.random.choice(target_neighbors,forward_burn_number, replace=False)
            
            #add more neighbors, with possible return edges
            for i in forward_burn_set:
                #outward edges
                G.add_edge(new_node_number,i)
                #now search the ones we just added
                list_of_targets.append(i)
        #end
    if(reverse_flag):
        G.reverse(copy=False)
    f = open(directory_path + file_name,'wb')
    pickle.dump(G,f)
    f.close()
    #nx.draw(G)
    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=20)
    nx.draw_networkx_edges(G, pos)
    pl.savefig('test.png')

    return G

def generate_gnr(n,p = 0.1,file_name='test_graph.p',directory_path='./cascade_data/'):
    G = nx.gnr_graph(n,p)
    f = open(directory_path + file_name,'wb')
    pickle.dump(G,f)
    f.close()
    #nx.draw(G)
    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=20)
    nx.draw_networkx_edges(G, pos)
    pl.savefig('test.png')
    return G

def generate_sc(n,alpha=.1,beta=.8, gamma=.1,file_name='test_graph.p',directory_path='./cascade_data/'):
    G = nx.scale_free_graph(n,alpha=alpha,beta=beta,gamma=gamma)
    #G = nx.scale_free_graph(n)
    f = open(directory_path+file_name,'wb')
    pickle.dump(G,f)
    f.close()
    #nx.draw(G)
    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=20)
    nx.draw_networkx_edges(G, pos)
    pl.savefig('test.png')
    return G

def generate_ba_graph(n,m=2,file_name='test_graph.p',directory_path='./cascade_data/'):
    G = nx.barabasi_albert_graph(n,m)
    f = open(directory_path+file_name,'wb')
    pickle.dump(G,f)
    f.close()
    #nx.draw(G)
    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=20)
    nx.draw_networkx_edges(G, pos)
    pl.savefig('test.png')
    return G

def generate_pc_graph(n,m=2, p=0.5,file_name='test_graph.p',directory_path='./cascade_data/'):
    G = nx.powerlaw_cluster_graph(n,m,p)
    f = open(directory_path+file_name,'wb')
    pickle.dump(G,f)
    f.close()
    #nx.draw(G)
    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=20)
    nx.draw_networkx_edges(G, pos)
    pl.savefig('test.png')
    return G

def generate_kron_graph(iterations = 1, seed = np.array([[0.962,0.107],[0.107,0.962]]), file_name='test_graph.p',directory_path='./cascade_data/'):
    from apgl.graph import *
    from apgl.generator.KroneckerGenerator import KroneckerGenerator
    
    initialGraph = SparseGraph(2)
    for i in range(2):
        for j in range(2):
            initialGraph[i,j] = seed[i,j]
            
    generator = KroneckerGenerator(initialGraph, iterations)
    graph = generator.generate()
    G = nx.DiGraph(graph.adjacencyMatrix())
    f = open(directory_path+file_name,'wb')
    pickle.dump(G,f)
    f.close()
    #nx.draw(G)
    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=20)
    nx.draw_networkx_edges(G, pos)
    pl.savefig('test.png')
    return G

if __name__ == '__main__':
    np.random.seed(1919)
    generate_pc_graph(50)
    #generate_ba_graph(100)
    #generate_sc(100)
    #generate_gnr(100)
    #generate_forest_fire_network(100,reverse_flag=False)