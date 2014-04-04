# -*- coding: utf-8 -*-
#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import glob
import os

def plot_measures(measures_files, k_max, name):
    '''Measures curves. 
    For each gyrus, plot measures values according to k clusters, 
    and the average. In each directory, a curves has been saved, 
    so for each gyrus.
    
    Parameters:
    ----------
        - measures_files (directory): contain .npy files for on gyrus
        - k_max (integer): the max number of cluster in the clustering
        - name (string): measure name
    '''
    for i in range(len(measures_files)):
        # finds all the pathnames matching a specified pattern
        inputs = glob.glob(str(measures_files[i]) + '/*' + name + '.npy')
        
        # usual colors
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
        
        # name of partition: needs to be more general
        partition = ['A', 'B', 'C']

        listof_inputs = []
        rand_index_fig = plt.figure()
        j = 0
        # plot measures values for each partition for one gyrus
        for infile in inputs:
            listof_inputs.append(np.loadtxt(infile))
            plt.plot(range(2, k_max + 1), np.loadtxt(infile), colors[j], 
                     label=("partition_" + partition[j]))
            plt.plot(range(2, k_max + 1), np.loadtxt(infile), 
                     colors[j] + 'o')
            plt.ylabel(name)
            plt.xlabel('Clusters (K)')
            # Mieux quand il y aura la table de correspondance
            gyrus = os.path.basename(str(measures_files[i]))
            rand_index_fig.suptitle(gyrus)
            plt.legend()
            j += 1
        
        # includes all arrays
        fig = np.vstack(listof_inputs)
        del listof_inputs
        
        # calculate the average for one gyrus on partitions
        mean_fig = np.mean(fig, axis=0)
        
        # plot the average
        plt.plot(range(2, k_max + 1), mean_fig, 'k--', label="average")
        plt.plot(range(2, k_max + 1), mean_fig, 'ko')
        plt.legend(loc='lower right')
        
        # save to svg image
        output_fig = os.path.join(str(measures_files[i]), name + ".svg")
        rand_index_fig.savefig(output_fig)