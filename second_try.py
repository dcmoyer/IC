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

def diff(a, b):
    b = set(b)
    return [aa for aa in a if aa not in b]

G = nx.Graph()
cascade_list = []
not_seen_list = []
cascade_adjacencies = []

def change_the_adj(adj_to_change, i,j, beta, epsilon, missing_p):
    time_evaluator = stats.expon.pdf
    delta = 0
    for c in adj_to_change:
        current_cascade = cascade_list[c]
        A = np.zeros(shape=(len(G.nodes()),len(G.nodes())))
        already_selected = [current_cascade[0][1]] #add index of source
        for index in range(len(current_cascade)):
            if(index == 0):
                continue
            for parent_index in range(index):
                for x in not_seen_list[c]:
                    #if the double jump occurs or half of the jump occurs and the otherhalf has just been added
                    if( ((current_cascade[parent_index][1],x) in G.edges() or ((j,i) == (current_cascade[parent_index][1],x))) \
                        and ((x,current_cascade[index][1]) in G.edges() or ((j,i) == (current_cascade[index][1],x)))):
                        # P(j,x)*P(x,i) = p \beta^2 * (\frac{1}{2} * \Delta_{j,i})
                        A[x,index] = missing_p * beta * beta\
                            * (time_evaluator((current_cascade[index][0] - current_cascade[parent_index][0])/2) ** 2)
                #end
                if( ((current_cascade[parent_index][1],current_cascade[index][1]) in G.edges()) or \
                    ((current_cascade[parent_index][1],current_cascade[index][1]) == (j,i))):
                    A[parent_index,index] = beta * time_evaluator(current_cascade[index][0] - current_cascade[parent_index][0])
                else:
                    A[parent_index,index] = epsilon * time_evaluator(current_cascade[index][0] - current_cascade[parent_index][0])
                
            #end
        #select parents
        #this sub-problem is NP-hard?
        for x in not_seen_list[c]:
            index = np.argmax(A[x,:])
            while(np.argmax(A[:,index]) != x and A[x,index] != 0):
                A[x,index] = 0
            if(A[x,index] == 0):
                continue
            A[x,:] = A[x,:] * np.zeros(shape=(1,len(A[x,:])))
        
        cascade_adjacencies[c] = A
    return

def find_next_link_with_adj_no_iter(input):
    j,i,beta, epsilon, missing_p = input
    time_evaluator = stats.expon.pdf
    delta = 0
    adjs_to_change = []
    for c in range(len(cascade_list)):
        #we split this into 4 cases
        if(j not in cascade_list[c].nodes()):
            #1) if the parent is not in the cascade
            if(i not in cascade_list[c].nodes()):
                #1a) if both the parent and the child are not in the cascade
                continue
            
            #1b) the child is in the cascade but the parent is not
            max_double_jump_weight = 0
            max_index = -1
            for x in G.predecessors(j):
                if(cascade_list[c].node[x]['time'] > cascade_list[c].node[i]['time']):
                    continue
                double_jump_weight = np.log(beta * beta * missing_p * \
                    ((time_evaluator(cascade_list[c].node[i]['time'] - cascade_list[c].node[x]['time'])/2)**2) ) - np.log(epsilon)
                if(double_jump_weight > max_double_jump_weight and double_jump_weight > np.max(cascade_adjacencies[c][:,i])):
                    max_index = x
                    max_double_jump_weight = double_jump_weight
            
            if(max_double_jump_weight > 0):
                delta = delta + max_double_jump_weight - np.max(cascade_adjacencies[c][:,i])
                adj_to_change.append(c)
        elif(i not in cascade_list[c].nodes()):
            #2) if the child is not in the cascade but the parent is
            max_double_jump_weight = 0
            max_index = -1
            for x in G.sucessors(i):
                if(cascade_list[c].node[x]['time'] < cascade_list[c].node[j]['time']):
                    continue
                double_jump_weight = beta * beta * missing_p * \
                    ((time_evaluator(cascade_list[c].node[x]['time'] - cascade_list[c].node[j]['time'])/2)**2) ) - np.log(epsilon)
            if(max_double_jump_weight > 0):
                delta = delta + max_double_jump_weight - np.max(cascade_adjacencies[c][:,i])
                adj_to_change.append(c)
        elif(cascade_list[c].node[i]['time'] < cascade_list[c].node[j]['time']):
            #3) both are in the cascade but the parent appears after the child
            continue
        elif(i in cascade_list[c].nodes() and j in cascade_list[c].nodes()):
            #4) j preceeds i, both in the cascades
            new_weight = np.log(beta * time_evaluator(cascade_list[c].node[i]['time'] - cascade_list[c].node[j]['time'])) - np.log(epsilon)
            old_weight = np.argmax(cascade_adjacencies[c][:,i])
            if(new_weight > old_weight)
                delta = delta + new_weight - old_weight 
                adj_to_change.append(c)
    return (delta, (j,i), adj_to_change)
    

def NetInf(k, graph_name, number_of_cascades, cascade_directory, missing_p, num_threads = 8, time_evaluator= stats.expon.pdf, beta = 0.5, epsilon = 0.0001):
    #
    #Initialization
    #
    global G
    G = nx.Graph()
    f = open(cascade_directory + graph_name + '.p','rb')
    G.add_nodes_from(pickle.load(f))
    f.close()
    
    global cascade_list
    cascade_list = []
    global not_seen_list
    not_seen_list = []
    #
    #Initial Trees(!)
    #
    global cascade_adjacencies
    
    for c in range(number_of_cascades):
        f = open(cascade_directory + graph_name + '/' +graph_name + '_' + str(c) + '.p', 'rb')
        current_cascade = pickle.load(f)
        f.close()
        
        cascade_list.append(current_cascade);
        cascade_adjacencies.append(np.zeros(size(len(G.nodes()))))
    
    for c in range(number_of_cascades):
        not_seen_list.append(diff(G.nodes(), [tuple[1] for tuple in cascade_list[c]]))
    
    
    #
    #Actual Algorithm
    #
    for i in range(k):
        print(i)
        #initial shifts
        #deltas = np.zeros(shape=(len(G.nodes())))
        
        max_delta = -np.inf
        max_pair = None
                
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
        block_size = len(active_set)//(num_threads) + 1
        #remainder = len(active_set) - (block_size * num_threads)
        
        #create oiik
        pool = Pool(num_threads)
        
        #tasks = [(pair[0],pair[0]+1, pair[1], pair[1]+1, beta, time_evaluator, tree_list) for pair in active_set]
        
        tasks = [(pair[0], pair[1], beta, epsilon, missing_p) for pair in active_set]
        output = pool.imap(find_next_link_BF_no_iter,tasks, block_size)
        pool.close()
        pool.join()
        for output_tuple in output:
            #output_tuple = output_result.get(timeout=1)
            if output_tuple[0] > max_delta:
                max_delta = output_tuple[0]
                max_pair = output_tuple[1]
                adj_to_change = adj_to_change
        print(max_delta)
        print(max_pair)
        
        change_the_adj(adj_to_change, i, j, beta, epsilon, missing_p)
        #end
        
        #ok so we found the max_delta
        #let us change the trees
        
        G.add_edge(max_pair[0],max_pair[1])
        return G
    #end
    
    f = open(cascade_directory + graph_name + '_NETINF_mult.p','wb')
    pickle.dump(G,f)
    f.close()
    
    return G
#end