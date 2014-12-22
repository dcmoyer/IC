
import matplotlib
matplotlib.use('Agg')

import numpy as np
import networkx as nx
import pickle
import matplotlib.pyplot as plt

def diff(a, b):
    b = set(b)
    return [aa for aa in a if aa not in b]

beta = 0.7
cascade_directory = './results/'
graphname = 'Core_Periphery_small'
graph_name = 'Core_Periphery_small_noise'
number_of_cascades = 100
p_values = [0.04]#[0.01,0.03,0.05,0.07,0.09]
missing_p_values = [0.2]#,0.5,1,2]

accuracy_array = np.zeros((len(p_values), len(missing_p_values)))
baseline_array = np.zeros((len(p_values),1))
for i,p in enumerate(p_values):
    
    ##BASELINE
    print('\nBASELINE ACCURACY:\n')
    print(p)
    f = open(cascade_directory + graph_name + str(p) + '_NETINF_mult.p','rb')
    G = pickle.load(f)
    f.close()
    
    f = open(cascade_directory + graph_name + str(p) + '.p','rb')
    E_objected = pickle.load(f).edges()
    f.close()
    
    baseline_array[i,0] = str(len(diff(G.edges(),E_objected))/float(len(E_objected)))
    
    #print('Type 2 error: ' + str(len(diff(G.edges(),E_objected))/float(len(E_objected))))
    print('Type 1 error: ' + str(len(diff(E_objected,G.edges()))/float(len(E_objected))))
    
    ##NEW##
    print('\nNew ACCURACY:\n')
    for j,missing_p in enumerate(missing_p_values):
        
        
        f = open(cascade_directory + graph_name + str(p) + '_newINF' + str(missing_p) + '_mult.p','rb')
        G = pickle.load(f)
        f.close()
        
        f = open(cascade_directory + graph_name + str(p) + '.p','rb')
        E_objected = pickle.load(f).edges()
        f.close()
        
        accuracy_array[i,j] = str(len(diff(G.edges(),E_objected))/float(len(E_objected)))
        print(missing_p)
        #print('Type 2 error: ' + str(len(diff(G.edges(),E_objected))/float(len(E_objected))))
        print('Type 1 error: ' + str(len(diff(E_objected,G.edges()))/float(len(E_objected))))

for i in range(len(missing_p_values)):
    fig = plt.figure()
    plt.plot(p_values, accuracy_array[:,i],'r')
    plt.plot(p_values, baseline_array[:,0],'b')
    plt.savefig( './plots/' + str(i) + '.png')