#!/usr/bin/python
#
#  IC
#
#

from generate_network import generate_ba_graph
from generate_cascades import generate_cascades
import netInf
import numpy as np
import networkx as nx
import pickle

beta = 0.7
cascade_directory = './cascade_data/'
graph_name = 'ba_100'
number_of_cascades = 1000

np.random.seed(1919)
generate_ba_graph(100,m=2,file_name='ba_100.p',directory_path=cascade_directory)
generate_cascades(number_of_cascades, spread_probability=beta , graph_file='ba_100.p',directory_path=cascade_directory)

f = open(cascade_directory + graph_name + '.p','rb')
k = len(pickle.load(f).edges())
f.close()

G = netInf.NetInf(k, graph_name, number_of_cascades, cascade_directory)


