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
*

Main dependencies:

Author:
"""


# ---------------------------Imports-------------------------------------------


# Soma module
from soma import aims
import six


# ---------------------------Functions-----------------------------------------


def merge_bundlegraph_and_roigraph(roi_graph,
                                   bundles_graph,
                                   bundles_to_ROIs_motion,
                                   nodes_names_mapping={}):
    """Merge Aims ROIs graph with Aims bundles graph.

    Parameters
    ----------
        roi_graph: roi_name (str),
                   roi_label (int),
                   roi_mesh (mesh),
                   aims_roi (bucket)
        bundles_graph:
        bundles_to_ROIs_motion: transformation matrix (dw to t1)
        nodes_names_mapping: (optional) dict
            bundles_graph to roi_graph nodes names correspondance map.
            Maps from bundles_graph "name" attribute to roi_graph "roi_label"
            attribute (which uses int values).
            Names not in this map are not translated and kept as is.

    Return
    -------
        roi_graph: add to roi_graph bundles as edges between nodes (=ROIs)
    """

    # associated ROIs
    roi_vertices = roi_graph.vertices().list()
    roi_vertices_labels = {}

    for i in six.moves.xrange(roi_vertices.__len__()):
        # if a given key is available, we keep the vertex
        # and roi_vertices[i]['roi_label'] != 0:
        if "roi_label" in roi_vertices[i]:
            roi_vertices_labels[
                int(roi_vertices[i]['roi_label'])] = roi_vertices[i]
    del roi_vertices

    # creating a trash graph
    # (What it is for, just to create roi_vertices_labels[0] ? )
    # -> for relations with the background node
    if 0 not in roi_vertices_labels:
        others_vertex = roi_graph.addVertex('roi')
        others_vertex['name'] = '0'
        roi_vertices_labels[0] = others_vertex.get()

    #other_name = nodes_names_mapping.get('others', 0)
    #others_vertex = None
    #if other_name:
        #others_vertex = roi_vertices_labels.get(other_name)
    #if others_vertex is None:
        #others_vertex = roi_graph.addVertex('roi')
        #roi_vertices_labels[0] = others_vertex.get()
        #others_vertex['name'] = '0'
    #others_vertex['name'] = 'others'

    # creating edges between labels/fibers
    for bundles in bundles_graph.vertices():
        names = bundles['name'].split('_')
        try:
            roi_1 = int(names[0])
            roi_2 = int(names[1])
        except ValueError:
            # a non-numeric name ('trash') has been used
            continue
        roi_1 = nodes_names_mapping.get(roi_1, roi_1)
        roi_2 = nodes_names_mapping.get(roi_2, roi_2)
        if roi_1 != roi_2:
            edge = roi_graph.addEdge(roi_vertices_labels[roi_1],
                                     roi_vertices_labels[roi_2],
                                     'roi_junction')
            fiber_mesh = bundles['aims_Tmtktri']
            aims.SurfaceManip.meshTransform(fiber_mesh.get(),
                                            bundles_to_ROIs_motion)
            aims.GraphManip.storeAims(roi_graph, edge.get(),
                                      'roi_mesh_junction', fiber_mesh)
