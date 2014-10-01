#!/usr/bin/python
#
#  IC
#
#

from generate_network import generate_ba_graph
from generate_cascades import generate_cascades
import netInf_mult
import netInf
import numpy as np
import networkx as nx
import pickle
import time
import first_try
import second_try
import noise

beta = 0.7
cascade_directory = './cascade_data/'
graph_name = 'ba_50'
number_of_cascades = 50

# np.random.seed(1919)
# generate_ba_graph(50,m=2,file_name=graph_name+'.p',directory_path=cascade_directory)
# generate_cascades(number_of_cascades, spread_probability=beta , graph_file='ba_50.p',directory_path=cascade_directory)
# noise.noise_the_data(graph_name,cascade_directory,number_of_cascades, 0.1)

f = open(cascade_directory + graph_name + '.p','rb')
k = len(pickle.load(f).edges())
f.close()

#start = time.time()
#G = netInf.NetInf(k, graph_name, number_of_cascades, cascade_directory,beta=beta)
#stop = time.time()
#print("Single Thread Time: %f" % (stop - start))

print("k = %d" % k)

start = time.time()
G = netInf_mult.NetInf(k, graph_name, number_of_cascades, cascade_directory,beta=beta,num_threads=6)
stop = time.time()

print("Multi Thread Time: %f" % (stop - start))

start = time.time()
second_try.NetInf(k, graph_name, number_of_cascades, cascade_directory, 0,beta=beta,num_threads=6)
stop = time.time()

print("New Algorithm Time: %f" % (stop - start))