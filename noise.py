#!/usr/bin/python
#
#noise.py
#
#Noise The Data

#import matplotlib
#matplotlib.use('Agg')

import networkx as nx
import numpy as np
import pickle
#import matplotlib.pyplot as pl
import csv
from scipy import stats

def noise_the_data(graph_name,cascade_directory,number_of_cascades , p= 0.01):
    
    for c in range(number_of_cascades):
        f = open(cascade_directory + graph_name + '/' +graph_name + '_' + str(c) + '.p', 'rb')
        current_cascade = pickle.load(f)
        f.close()
        
        new_cascade = []
        for x in current_cascade:
            if(np.random.random() > p):
                new_cascade.append(x)
        
        f = open(cascade_directory + graph_name + '_noise' + str(p) + '/' + graph_name + '_noise' + str(p)+ '_' + str(c) + '.p', 'wb')
        pickle.dump(new_cascade,f)
        f.close()
    
    return
    

if __name__ == '__main__':
    noise_the_data('ba_50','./cascade_data/',500)