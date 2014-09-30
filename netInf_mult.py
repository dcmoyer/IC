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
from multiprocessing import Pool
import time

def diff(a, b):
    b = set(b)
    return [aa for aa in a if aa not in b]

def initialize_DAG(cascade,epsilon):
    T = nx.DiGraph()
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

def find_next_link(j_start,j_end,i_start,i_end,beta,time_evaluator,tree_list):
    # if(i==j):
    # return (0,[])
    #don't add edges already in the graph
    # if((j,i) in G.edges()):
    # return (0,[])
    listing = []
    for j in range(j_start,j_end):
        for i in range(i_start,i_end):
            delta = 0
            trees_to_change_temp = []
            #deltas[i,j] = 0 #implicit
            for c in range(len(tree_list)):
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
                    
                    #calculate the new weight
                    new_weight = np.log(beta * time_evaluator(tree_list[c].node[i]['time'] - tree_list[c].node[i]['time'])) - np.log(epsilon)
                    #if the new weight is bigger than the old
                    old_weight = tree_list[c][tree_list[c].node[i]['parent']][i]['weight']
                    if(new_weight >= old_weight):
                        
                        delta = delta + new_weight - old_weight
                        trees_to_change_temp.append((c,new_weight))
                        
                    else:
                        continue
                else: #neither active
                    continue
            #end
            listing.append((delta, (j,i),trees_to_change_temp))
        #end
    #end
    return listing

def find_next_link_no_iter(input):
    j,i,beta, epsilon,tree_list = input
    #print("%d, %d being evaluated" % (j,i))
    time_evaluator = stats.expon.pdf
    # if(i==j):
    # return (0,[])
    #don't add edges already in the graph
    # if((j,i) in G.edges()):
    # return (0,[])
    
    delta = 0
    trees_to_change_temp = []
    for c in range(len(tree_list)):
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
            
            #calculate the new weight
            new_weight = np.log(beta * time_evaluator(tree_list[c].node[i]['time'] - tree_list[c].node[j]['time'])) - np.log(epsilon)
            #if the new weight is bigger than the old
            old_weight = tree_list[c][tree_list[c].node[i]['parent']][i]['weight']
            if(new_weight >= old_weight):
                
                delta = delta + new_weight - old_weight
                trees_to_change_temp.append((c,new_weight))
                
            else:
                continue
        else: #neither active
            continue
    #end
    return (delta, (j,i),trees_to_change_temp)
    
def NetInf(k, graph_name, number_of_cascades, cascade_directory, num_threads = 8, time_evaluator= stats.expon.pdf, beta = 0.5, epsilon = 0.0001):
    #
    #Initialization
    #
    G = nx.DiGraph()
    f = open(cascade_directory + graph_name + '.p','rb')
    G.add_nodes_from(pickle.load(f))
    f.close()
    
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
    overall_start = time.time()
    for i in range(k):
        print(i)
        start = time.time()
        #initial shifts
        #deltas = np.zeros(shape=(len(G.nodes())))
        
        max_delta = -np.inf
        max_pair = None
        
        trees_to_change_max = []
        
        active_set = []
        
        for j,i in [(x,y) for x in range(len(G.nodes())) for y in range(len(G.nodes()))]:
            #no self loops
            if(i==j):
                continue
            #don't add edges already in the graph
            if((j,i) in G.edges()):
                continue
            active_set.append((j,i))
        
        #determine block size
        block_size = len(active_set)//(num_threads)
        #remainder = len(active_set) - (block_size * num_threads)
        
        #create oiik
        pool = Pool(num_threads)
        
        #tasks = [(pair[0],pair[0]+1, pair[1], pair[1]+1, beta, time_evaluator, tree_list) for pair in active_set]
        
        tasks = [(pair[0], pair[1], beta, epsilon, tree_list) for pair in active_set]
        output = pool.imap(find_next_link_no_iter,tasks, block_size)
        pool.close()
        pool.join()
        
        for output_tuple in output:
            #output_tuple = output_result.get(timeout=1)
            if output_tuple[0] > max_delta:
                max_delta = output_tuple[0]
                trees_to_change_max = output_tuple[2]
                max_pair = output_tuple[1]
        
        print(max_delta)
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
        stop = time.time()
        print("Iteration time = %f" % (stop - start))
        print("Overall time = %f" % (stop - overall_start))
    #end
    
    f = open(cascade_directory + graph_name + '_NETINF_mult.p','wb')
    pickle.dump(G,f)
    f.close()
    
    return G
#end