#!/usr/bin/python
#
#

import matplotlib
matplotlib.use('Agg')

import networkx as nx
import numpy as np
import pickle
import matplotlib.pyplot as pl
import Queue


def diff(a, b):
    b = set(b)
    return [aa for aa in a if aa not in b]

def generate_cascades(n, spread_probability ,time_generator=np.random.exponential,\
 graph_file='test_graph.p',directory_path='./cascade_data/'):
    f = open(directory_path + graph_file,'rb')
    G = pickle.load(f)
    
    for cascade_index in range(n):
        
        #choose a source
        s = int(np.random.choice(G.nodes(),1))
        
        already_seen = [s]
        
        Q = Queue.PriorityQueue()
        Q.put((0,s))
        
        current_cascade = []
        
        while(True):
            if Q.empty():
                break
            next_node = Q.get()
            node_id = next_node[1]
            current_time = next_node[0]
            current_cascade.append((next_node[0],next_node[1]))
            
            #decide on neighbors to spread to
            unseen_neighbors = diff(G.neighbors(node_id),already_seen)
            for i in unseen_neighbors:
                if np.random.uniform() < spread_probability:
                    Q.put((current_time + time_generator(),i))
                    already_seen.append(i)
                #end
            #end
            
        #end
        file_output = open(directory_path + graph_file[:-2] + '/' + graph_file[:-2] +'_' + str(cascade_index) + '.p','wb')
        pickle.dump(current_cascade,file_output)
        file_output.close()
    #end
    file_output = open(directory_path + graph_file[:-2] + '/' + graph_file[:-2] + '_node_list.p','wb')
    pickle.dump(G.nodes(),file_output)
    file_output.close()
    f.close()

if __name__ == "__main__":
    generate_cascades(1, spread_probability=0.5)