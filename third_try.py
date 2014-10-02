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

G = nx.DiGraph()
cascade_list = []
not_seen_list = []
cascade_adjacencies = []
timings = []
LOG_EPSILON = 1

def change_the_adj(input):
    c, i,j, beta, epsilon, missing_p = input
    time_evaluator = stats.expon.pdf
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
                    A[x,current_cascade[index][1]] = missing_p * beta * beta\
                        * (time_evaluator((current_cascade[index][0] - current_cascade[parent_index][0])/2) ** 2)
            #end
            if( ((current_cascade[parent_index][1],current_cascade[index][1]) in G.edges()) or \
                ((current_cascade[parent_index][1],current_cascade[index][1]) == (j,i))):
                A[current_cascade[parent_index][1],current_cascade[index][1]] = beta * time_evaluator(current_cascade[index][0] - current_cascade[parent_index][0])
            else:
                A[current_cascade[parent_index][1],current_cascade[index][1]] = epsilon * time_evaluator(current_cascade[index][0] - current_cascade[parent_index][0])
        #end
    #select parents
    #this sub-problem is NP-hard?
    '''for x in not_seen_list[c]:
        index = np.argmax(A[x,:])
        while(np.argmax(A[:,index]) != x and A[x,index] != 0):
            A[x,index] = 0
        if(A[x,index] == 0):
            continue
        A[x,:] = A[x,:] * np.zeros(shape=(1,len(A[x,:])))
    '''
    return A

def find_next_link_with_adj_no_iter(input):
    j,i,beta, epsilon, missing_p = input
    time_evaluator = stats.expon.pdf
    delta = 0
    adj_to_change = []
    for c in range(len(cascade_list)):
        #we split this into 4 cases
        if(j in not_seen_list[c]):
            #1) if the parent is not in the cascade
            if(i in not_seen_list[c]):
                #1a) if both the parent and the child are not in the cascade
                continue
            
            #1b) the child is in the cascade but the parent is not
            max_double_jump_weight = 0
            max_index = -1
            current_max = np.log(np.max(cascade_adjacencies[c][:,i])) - LOG_EPSILON
            for x in G.predecessors(j):
                if(x in not_seen_list[c]):
                    continue
                if(timings[c][x] > timings[c][i]):
                    continue
                double_jump_weight = np.log(beta * beta * missing_p * \
                    ((time_evaluator(timings[c][i] - timings[c][x])/2)**2) ) - LOG_EPSILON
                if(double_jump_weight > max_double_jump_weight and double_jump_weight > current_max):
                    max_index = x
                    max_double_jump_weight = double_jump_weight
            if(max_double_jump_weight > 0):
                delta = delta + max_double_jump_weight - current_max
                adj_to_change.append(c)
        elif(i in not_seen_list[c]):
            #2) if the child is not in the cascade but the parent is
            max_double_jump_weight = 0
            max_index = -1
            for x in G.successors(i):
                if(x in not_seen_list[c]):
                    continue
                if(timings[c][x] < timings[c][j]):
                    continue
                current_max = np.log(np.max(cascade_adjacencies[c][:,x])) - LOG_EPSILON
                double_jump_weight = (beta * beta * missing_p * \
                    ((time_evaluator(timings[c][x] - timings[c][j])/2)**2))  - LOG_EPSILON
                if(double_jump_weight > max_double_jump_weight and double_jump_weight > current_max):
                    max_index = x
                    max_double_jump_weight = double_jump_weight
            if(max_double_jump_weight > 0):
                delta = delta + max_double_jump_weight - current_max
                adj_to_change.append(c)
        elif(timings[c][i] < timings[c][j]):
            #3) both are in the cascade but the parent appears after the child
            continue
        elif(i not in not_seen_list[c] and j not in not_seen_list[c]):
            #4) j preceeds i, both in the cascades
            new_weight = np.log(beta * time_evaluator(timings[c][i] - timings[c][j])) - LOG_EPSILON
            old_weight = np.log(np.max(cascade_adjacencies[c][:,i])) - LOG_EPSILON
            if(new_weight > old_weight):
                delta = delta + new_weight - old_weight 
                adj_to_change.append(c)
    return (delta, (j,i), adj_to_change)
    

def NetInf(k, graph_name, number_of_cascades, cascade_directory, missing_p, verbose=False, num_threads = 8, time_evaluator= stats.expon.pdf, beta = 0.5, epsilon = 0.0001):
    #
    #Initialization
    #
    global G
    G = nx.DiGraph()
    f = open(cascade_directory + graph_name + '.p','rb')
    G.add_nodes_from(pickle.load(f))
    f.close()
    global LOG_EPSILON
    LOG_EPSILON = np.log(epsilon)
    global cascade_list
    cascade_list = []
    global not_seen_list
    not_seen_list = []
    #
    #Initial Trees(!)
    #
    global cascade_adjacencies
    global timings
    
    overall_start = time.time()
    
    for c in range(number_of_cascades):
        f = open(cascade_directory + graph_name + '/' +graph_name + '_' + str(c) + '.p', 'rb')
        current_cascade = pickle.load(f)
        f.close()
        if current_cascade == []:
            number_of_cascades = number_of_cascades - 1
            continue
        cascade_list.append(current_cascade);
        cascade_adjacencies.append(np.zeros(shape=(len(G.nodes()),len(G.nodes()))))
    
    for c in range(number_of_cascades):
        not_seen_list.append(diff(G.nodes(), [tuple[1] for tuple in cascade_list[c]]))
    
    for c in range(number_of_cascades):
        timings.append({node_id : time for time,node_id in cascade_list[c]})
    #
    #Actual Algorithm
    #
    block_size = number_of_cascades//(num_threads) + 1
    pool = Pool(num_threads)
    tasks = [(c, -1,-1, beta, epsilon, missing_p) for c in range(number_of_cascades)]
    output = pool.imap(change_the_adj, tasks,block_size)
    pool.close()
    pool.join()
    for adj in range(number_of_cascades):
        cascade_adjacencies[adj] = output.next()
    #change_the_adj(range(number_of_cascades), -1, -1, beta, epsilon, missing_p)
    for step in range(k):
        
        print(step)
        start = time.time()
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
            if((j,i) in G.edges() or (i,j) in G.edges()):
                continue
            active_set.append((j,i))
        
        #determine block size
        block_size = len(active_set)//(num_threads) + 1
        #remainder = len(active_set) - (block_size * num_threads)
        
        #create oiik
        pool = Pool(num_threads)
        
        #tasks = [(pair[0],pair[0]+1, pair[1], pair[1]+1, beta, time_evaluator, tree_list) for pair in active_set]
        
        tasks = [(pair[0], pair[1], beta, epsilon, missing_p) for pair in active_set]
        output = pool.imap(find_next_link_with_adj_no_iter,tasks, block_size)
        #output = map(find_next_link_with_adj_no_iter, tasks)
        pool.close()
        pool.join()
        for output_tuple in output:
            #output_tuple = output_result.get(timeout=1)
            if output_tuple[0] > max_delta:
                max_delta = output_tuple[0]
                max_pair = output_tuple[1]
                adj_to_change = output_tuple[2]
        
        block_size = len(adj_to_change)//(num_threads) + 1
        pool = Pool(num_threads)
        tasks = [(c, max_pair[1], max_pair[0], beta, epsilon, missing_p) for c in adj_to_change]
        output = pool.imap(change_the_adj, tasks,block_size)
        #output = map(change_the_adj, tasks)
        pool.close()
        pool.join()
        
        for adj in range(len(adj_to_change)):
            cascade_adjacencies[adj_to_change[adj]] = output.next()
        
        #change_the_adj(adj_to_change, max_pair[1], max_pair[0], beta, epsilon, missing_p)
        #end
        
        #ok so we found the max_delta
        #let us change the trees
        
        G.add_edge(max_pair[0],max_pair[1])
        stop = time.time()
        if verbose:
            print(max_delta)
            print(max_pair)
            print("Iteration time = %f" % (stop - start))
            print("Overall time = %f" % (stop - overall_start))
    #end
    if not verbose:
        print("Overall time = %f" % (stop - overall_start))
    f = open(cascade_directory + graph_name + '_newINF_mult.p','wb')
    pickle.dump(G,f)
    f.close()
    
    return G
#end