#!/usr/bin/env python
###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################
import numpy as np
import random


def roi_graph(mesh, texture):

    # create a np array from gyri segmentation (aims to np)
    tex = texture[0].arraydata()

    triangles = np.asarray([tri.arraydata() for tri in mesh.polygon()])

    graph = {}
    for vertex in range(tex.size):
        vertex_neighbors = []

        # return the index of the triangles associated with a ref vertex
        idx_list = np.where(triangles == vertex)[0]

        # create a list with all the vertices connected to the ref vertex
        for idx in idx_list:
            vertex_neighbors.extend(triangles[idx])

        # remove multiple occurences of a vertex and the ref vertex
        vertex_neighbors = np.unique(vertex_neighbors)
        vertex_neighbors = vertex_neighbors[vertex_neighbors != vertex]

        # create the label array of the neighbors of a ref vertex
        labels = []
        for neighbor in vertex_neighbors:
            labels.append(tex[neighbor])
        labels = np.unique(labels)
        labels = labels[labels != tex[vertex]].tolist()

        # store the regions neighboring the region of a vertex
        if tex[vertex] in graph.keys():
            graph[tex[vertex]].extend(labels)
            simplified_list = np.unique(graph[tex[vertex]]).tolist()
            graph[tex[vertex]] = simplified_list
        else:
            graph[tex[vertex]] = labels

    return graph


def saturation(node, graph, node_colors):
    sat = 0
    for label in graph[node]:
        if node_colors[label] != 0:
            sat += 1
    return sat


def max_saturation(graph, node_colors):
    sorted_nodes = sorted(graph, key=lambda node: -len(graph[node]))
    max = -1
    max_sat = -1
    for node in sorted_nodes:
        node_sat = saturation(node, graph, node_colors)
        if node_sat > max_sat:
            max = node
            max_sat = node_sat
    return max


def available_color(node, graph, node_colors, nb_colors="minimal"):
    colors = [1, 2, 3, 4, 5, 6, 7, 8]
    if nb_colors != "minimal":
        colors = colors[:int(nb_colors)]

    for label in graph[node]:
        if node_colors[label] in colors:
            colors.remove(node_colors[label])

    if nb_colors == "minimal":
        return colors[0]
    else:
        random_color = random.randint(0, len(colors)-1)
        return colors[random_color]


def dsatur(graph, nb_colors="minimal"):
    node_colors = dict.fromkeys(graph.keys(), 0)
    while 0 in node_colors.values():
        node = max_saturation(graph, node_colors)
        color = available_color(node, graph, node_colors, nb_colors)
        node_colors[node] = color
        del graph[node]

    return node_colors


def get_colormap_colors(mesh, texture, nb_colors="minimal", default_nodes=[]):

    # list of colors : Default, Red, Green, Blue, Yellow, Purple, light Red,
    # light Green, Grey
    colors = {0: [255, 255, 255, 255],
              1: [203, 0, 0, 255],
              2: [64, 173, 38, 255],
              3: [0, 0, 142, 255],
              4: [255, 217, 0, 255],
              5: [186, 85, 211, 255],
              6: [203, 80, 90, 255],
              7: [60, 100, 20, 255],
              8: [100, 100, 100, 255]}

    # Create the graph of the parcellation on the mesh
    graph = roi_graph(mesh, texture)

    # Get color of labels with DSATUR algorithm
    node_colors = dsatur(graph, nb_colors)

    sorted_nodes = sorted(node_colors.keys())
    last_node = sorted_nodes[0]
    RGBA_colors = []
    for node in sorted_nodes:

        RGBA_colors.extend(colors[0]*(node - last_node - 1))

        if node in default_nodes:
            RGBA_colors.extend(colors[0])
        else:
            RGBA_colors.extend(colors[node_colors[node]])
        last_node = node

    return RGBA_colors
