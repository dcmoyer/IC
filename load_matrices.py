#!/usr/bin/python
#
#  IC
#
#


import matplotlib
matplotlib.use('Agg')
import numpy as np
from scipy import sparse
import networkx as nx
import pickle
import matplotlib.pyplot as pl


def load_matrices(filename,directory_path='./cascade_data/'):
    f = open('xinrans_graphs/' + filename + '.txt','r')
    
    rows = []
    cols = []
    vals = []
    #tupes = []
    for line in f:
        line = line.rstrip()
        if line == '':
            continue
        row,col = map(int,line.split(','))
        
        #tupes.append((1,(row,col)))
        rows.append(row)
        cols.append(col)
        vals.append(1)
        
    G = nx.DiGraph(sparse.coo_matrix((vals,(rows,cols))))
    f = open(directory_path+filename + '.p','wb')
    pickle.dump(G,f)
    f.close()
    #nx.draw(G)
    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=20)
    nx.draw_networkx_edges(G, pos)
    pl.savefig('test.png')
    return G

if __name__ == '__main__':
    load_matrices('Core_Periphery_small')