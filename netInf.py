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
import csv
from scipy import stats

def diff(a, b):
    b = set(b)
    return [aa for aa in a if aa not in b]

def initialize_DAG(cascade,epsilon):
    T = nx.Graph()
    node_list = [c[1] for c in cascade]
    parent_id = -1
    for pair in cascade:
        T.add_node(pair[1],time=pair[0], parent=parent_id)
        parent_id = pair[1]
    for node_index in range(len(node_list)):
        if node_index == 0:
            continue
        T.add_edge(node_list[node_index-1],node_list[node_index],weight=0)
    return T

#TODO Make brute force
def bf_DAG(cascade):
    return

def NetInf(k, graph_name, number_of_cascades, cascade_directory, time_evaluator= stats.expon.pdf, beta = 0.5, epsilon = 0.0001):
    #
    #Initialization
    #
    G = nx.Graph()
    f = open(cascade_directory + graph_name + '.p','rb')
    G.add_nodes_from(pickle.load(f))
    f.close()
    
    # epsilon_edge_list = []
    # for i in range(len(G.nodes()):
    #     for j in range(len(G.nodes()):
    #         if( i == j ):
    #             continue
    #         epsilon_edge_list.append((i,j,epsilon))
    # G.add_weighted_edges_from(epsilon_edge_list)
    
    tree_list = []
    cascade_list = []
    
    #
    #Initial Trees(!)
    #
    for c in range(number_of_cascades):
        f = open(cascade_directory + graph_name + '/' +graph_name + '_' + str(c) + '.p', 'rb')
        current_cascade = pickle.load(f)
        f.close()
        
        cascade_list.append(current_cascade);
        
        tree_list.append(initialize_DAG(current_cascade,epsilon))
        
    #
    #Actual Algorithm
    #
    for i in range(k):
        print(i)
        #initial shifts
        #deltas = np.zeros(shape=(len(G.nodes())))
        
        max_delta = -np.inf
        max_pair = None
        
        trees_to_change_max = []
        for j,i in [(x,y) for x in range(len(G.nodes())) for y in range(len(G.nodes()))]:
            delta = 0
            #no self loops
            
            if(i==j):
                continue
            #don't add edges already in the graph
            if((j,i) in G.edges()):
                continue
            
            trees_to_change_temp = []
            #deltas[i,j] = 0 #implicit
            for c in range(number_of_cascades):
                #if i preceeds j in this cascade
                if(j not in tree_list[c].nodes()):
                    continue
                    #do nothing
                elif(i not in tree_list[c].nodes()): # i inactive
                    continue
                    #delta = delta \
                    #    + np.log(1 - beta) \
                    #    - np.log(1 - epsilon)
                    #we're using log of the edge weights
                    #why?
                    #   Because this makes us able to sum the probability,
                    #   instead of multiplying things together.
                    #we don't need trees_to_change_temp.append(c)
                    #because it's not active!
                elif((tree_list[c].node[i])['time'] < (tree_list[c].node[j])['time']):
                    continue
                elif (i in tree_list[c].nodes() and j in tree_list[c].nodes()):
                    
                    the_active_tree = tree_list[c]
                    
                    #calculate the new weight
                    new_weight = np.log(beta * time_evaluator(the_active_tree.node[i]['time'] - the_active_tree.node[i]['time'])) - np.log(epsilon)
                    #if the new weight is bigger than the old
                    old_weight = the_active_tree[the_active_tree.node[i]['parent']][i]['weight']
                    if(new_weight >= old_weight):
                        
                        delta = delta + new_weight - old_weight
                        trees_to_change_temp.append((c,new_weight))
                        
                    else:
                        continue
                        #otherwise the old is bigger than the new
                        #so we do nothing
                        #if the edge probabilites are different
                        #we may need to change this
                    #end
                #elif( i in tree_list[c].nodes()): # j inactive
                #    continue
                #    #nothing happens
                #    #as the weight is zero
                #    #so the new node is the same as the old node
                else: #neither active
                    continue
                
            #end
            
            if delta > max_delta:
                max_delta = delta
                trees_to_change_max = trees_to_change_temp
                max_pair = (j,i)
                print(delta)
                print(max_pair)
            #end
        
        #ok so we found the max_delta
        #let us change the trees
        
        for c,new_weight in trees_to_change_max:
            tree_list[c].remove_edge(tree_list[c].node[max_pair[1]]['parent'],max_pair[1])
            tree_list[c].add_edge(max_pair[0],max_pair[1], weight=new_weight)
            tree_list[c].node[max_pair[1]]['parent'] = max_pair[0]
        #end
        G.add_edge(max_pair[0],max_pair[1])
        
    #end
    
    f = open(cascade_directory + graph_name + '_NETINF.p','wb')
    pickle.dump(G,f)
    f.close()
    
    return G
#end