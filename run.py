#!/usr/bin/python
#
#  IC
#
#

#from generate_network import generate_ba_graph
from generate_cascades import generate_cascades
import netInf_mult
#import netInf
import numpy as np
import networkx as nx
import pickle
import time
#import first_try
#import second_try
import third_try
import noise

beta = 0.7
cascade_directory = './cascade_data/'
graphname = 'Core_Periphery_small'
graph_name = 'Core_Periphery_small_noise'
number_of_cascades = 100
p_values = [0.01,0.03,0.05,0.07,0.09,0.1]
missing_p_values = [0.2,0.5,1,2]
for p in p_values:
    for missing_p in missing_p_values:
        graph_name = 'Core_Periphery_small_noise' + str(p)
        
        # np.random.seed(1919)
        #generate_ba_graph(50,m=2,file_name=graph_name+'.p',directory_path=cascade_directory)
        generate_cascades(number_of_cascades, spread_probability=beta , graph_file=graphname + '.p',directory_path=cascade_directory)
        noise.noise_the_data(graphname,cascade_directory,number_of_cascades, p)
        #exit()
        
        f = open(cascade_directory + graph_name + '.p','rb')
        k = len(pickle.load(f).edges())
        f.close()
        
        #start = time.time()
        #G = netInf.NetInf(k, graph_name, number_of_cascades, cascade_directory,beta=beta)
        #stop = time.time()
        #print("Single Thread Time: %f" % (stop - start))
        print('graph_name')
        print("k = %d" % k)
        
        start = time.time()
        G = netInf_mult.NetInf(k, graph_name, number_of_cascades, cascade_directory,beta=beta,num_threads=8)
        stop = time.time()
        
        print("Multi Thread Time: %f" % (stop - start))
        
        start = time.time()
        third_try.NetInf(k, graph_name, number_of_cascades, cascade_directory, missing_p ,beta=beta,num_threads=8,verbose=False)
        stop = time.time()
        
        print("New Algorithm Time: %f" % (stop - start))
        if(missing_p == 0.2):
            start = time.time()
            third_try.NetInf(k, graph_name, number_of_cascades, cascade_directory, p ,beta=beta,num_threads=8,verbose=False)
            stop = time.time()
            print("New Algorithm Time: %f" % (stop - start))
    #end
#end