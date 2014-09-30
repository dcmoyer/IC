#errorcheck


from generate_network import generate_ba_graph
from generate_cascades import generate_cascades
import netInf
import numpy as np
import networkx as nx
import pickle

beta = 0.7
cascade_directory = './cascade_data/'
graph_name = 'ba_100'
number_of_cascades = 100

def diff(a, b):
    b = set(b)
    return [aa for aa in a if aa not in b]


f = open(cascade_directory + graph_name + '_NETINF_mult.p','rb')
G = pickle.load(f)
f.close()

f = open(cascade_directory + graph_name + '.p','rb')
E_objected = pickle.load(f).edges()
f.close()
print(len(diff(G.edges(), E_objected)) + len(diff(E_objected,G.edges())))

print(diff(G.edges(),E_objected))
print(diff(E_objected,G.edges()))

print('Type 2 error: ' + str(len(diff(G.edges(),E_objected))/float(len(E_objected))))
print('Type 1 error: ' + str(len(diff(E_objected,G.edges()))/float(len(E_objected))))
#print(E_objected)
#print(G.edges())