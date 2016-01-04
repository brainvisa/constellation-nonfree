#!/usr/bin/env python
###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################

"""
This script does the following:
*create pie chart
*

Main dependencies:

Author: Sandrine Lefranc, 2015
"""

#----------------------------Imports-------------------------------------------


#system import
import os
import numpy
from operator import itemgetter, add

# matplotlib import
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# soma import
from soma import aims

# constel import
from constel.lib.texturetools import identify_patch_number


#----------------------------Functions-----------------------------------------


def create_pie_chart(fname):
    """
    """
    # identify the number of the region from a filename
    region = identify_patch_number(fname)

    # read tracts files
    bundle = aims.read(fname)

    # the total number of tracts
    tot_nbtracts = bundle["curves_count"]

    # list of combinations of relations between regions R to be quantified
    # each relation r_r is associated with the first index to start accounting
    # for the tracts
    ltracts = bundle["bundles"]

    #
    relations = []
    nbtracts = []
    for i in range(len(ltracts)):
        if i % 2 == 0:
            relations.append(ltracts[i])
        else:
            nbtracts.append(ltracts[i])
    nb = []
    for i in range(len(nbtracts)):
        if nbtracts[i] == 0:
            n = nbtracts[i + 1]
        elif i == len(nbtracts) - 1:
            n = tot_nbtracts - nbtracts[i]
        else:
            n = nbtracts[i + 1] - nbtracts[i]
        nb.append(n)

    dico = {}
    for i in range(len(relations)):
        dico[relations[i]] = nb[i]

    return dico, tot_nbtracts, region


def remove_border(axes=None, top=False, right=False, left=False, bottom=True):
        """
        Minimize chartjunk by stripping out unnecesasry plot borders and axis
        ticks.The top/right/left/bottom keywords toggle whether the
        corresponding plot border is drawn
        """
        ax = axes or plt.gca()
        ax.spines['top'].set_visible(top)
        ax.spines['right'].set_visible(right)
        ax.spines['left'].set_visible(left)
        ax.spines['bottom'].set_visible(bottom)

        #turn off all ticks
        ax.yaxis.set_ticks_position('none')
        ax.xaxis.set_ticks_position('none')

        #now re-enable visibles
        if top:
            ax.xaxis.tick_top()
        if bottom:
            ax.xaxis.tick_bottom()
        if left:
            ax.yaxis.tick_left()
        if right:
            ax.yaxis.tick_right()


#----------------------------Main program--------------------------------------


if __name__ == "__main__":

    tfile = ("/neurospin/archi-public/Users/lefranc/bipolar/bipolarcontrols"
             "/subjects/ab100404/diffusion/default_acquisition"
             "/default_analysis/default_tracking_session/"
             "connectivity_parcellation/avg/fs27controls1/G27/filteredTracts"
             "/ab100404_fs27controls1_G27_fibersNearCortex_30to500mm.bundles")

    outdir = "/volatile/sandrine/test_visu"
    pdffile = os.path.join(outdir, "test_visu.pdf")
    tri = True

    dico, tot_nbtracts, patch = create_pie_chart(tfile)

    pp = PdfPages(pdffile)

    for element in dico.iteritems():
        print element[0][:2]
        if element[0][:2] != patch:
            pass

    if tri:
        d = sorted(dico.iteritems(), reverse=False, key=itemgetter(0))
        d = reduce(add, d)
        labels = []
        nb = []
        for i in range(len(d)):
            if i % 2 == 0:
                labels.append(d[i])
            else:
                nb.append(d[i])
        ypos = numpy.arange(len(labels))
        plt.title('Fiber tracts distribution', fontsize=20)
        plt.barh(ypos, nb, color="r", height=0.8)

        #add the numbers to the side of each bar
        for p, c, ch in zip(ypos, labels, nb):
            plt.annotate(
                str(int(ch)), xy=(ch + 5, p + 0.4), va='center', fontsize=6)

        #cutomize ticks
        plt.yticks(ypos + 0.5, labels, fontsize=6)
        xt = plt.xticks()[0]
        plt.xticks(xt, [' '] * len(xt))
        plt.xlabel(
            "Number of fiber tracts connected from region " + str(patch) +
            " \n total: " + str(int(tot_nbtracts)) + " fiber tracts",
            fontsize=14)
        plt.ylabel("Relation between patches", fontsize=14)

        #minimize chartjunk
        remove_border(left=False, bottom=False)
        plt.grid(axis='x', color='white', linestyle='-')

        #set plot limits
        plt.ylim(ypos.max() + 1, ypos.min() - 1)
        plt.xlim(0, max(nb) + 100)
        del labels, nb, ypos, d

    else:
        ypos = numpy.arange(len(dico))
        plt.title('Fiber tracts distribution', fontsize=18)
        plt.barh(ypos, dico.values(), color="r", height=0.8)

        #add the numbers to the side of each bar
        for p, c, ch in zip(ypos, dico.keys(), dico.values()):
            plt.annotate(
                str(int(ch)), xy=(ch + 10, p + 0.4), va='center', fontsize=6)

        #cutomize ticks
        plt.yticks(ypos + 0.5, dico.keys(), fontsize=6)
        xt = plt.xticks()[0]
        plt.xticks(xt, [' '] * len(xt))
        plt.xlabel(
            "Number of fiber tracts connected from region " + str(patch) +
            " \n total: " + str(int(tot_nbtracts)) + " fiber tracts",
            fontsize=14)
        plt.ylabel("Relation between patches", fontsize=14)

        #minimize chartjunk
        remove_border(left=False, bottom=False)
        plt.grid(axis='x', color='white', linestyle='-')

        #set plot limits
        plt.ylim(ypos.max() + 1, ypos.min() - 1)
        plt.xlim(0, max(dico.values()) + 100)

    # Done with the page
    pp.savefig()
    # Write the PDF document to the disk
    pp.close()
