# -*- coding: utf-8 -*-

# BrainVisa module
from soma import aims

def mergeBundlesGraphAndROIsGraph(roi_graph, bundles_graph, bundles_to_ROIs_motion):
    '''Merge Aims ROIs graph with Aims bundles graph.
    
    Parameters:
        roi_graph: roi_name (str), 
                   roi_label (int), 
                   roi_mesh (mesh), 
                   aims_roi (bucket)
        bundles_graph:
        bundles_to_ROIs_motion: transformation matrix (dw to t1)
    
    Result:
        roi_graph: add to roi_graph bundles as edges between nodes (=ROIs)
    '''
    
    # associated ROIs
    roi_vertices = roi_graph.vertices().list()
    roi_vertices_labels = {}
    
    for i in xrange(roi_vertices.__len__()):
        # if a given key is available, we keep the vertex
        if roi_vertices[i].has_key('roi_label') and roi_vertices[i]['roi_label'] != 0:
            roi_vertices_labels[int(roi_vertices[i]['roi_label'])] = roi_vertices[i]
    del roi_vertices
    
    # creating a trash graph 
    # (What it is for, just to create roi_vertices_labels[0] ? )
    others_vertex = roi_graph.addVertex('roi')
    others_vertex['name'] = 'others'
    roi_vertices_labels[0] = others_vertex.get()
    
    # creating edges between labels/fibers
    for bundles in bundles_graph.vertices():
        names = bundles['name'].split('_')
        roi_1 = int(names[0])
        roi_2 = int(names[1])
        if roi_1 != roi_2:
            edge = roi_graph.addEdge(roi_vertices_labels[roi_1], 
                                     roi_vertices_labels[roi_2], 
                                     'roi_junction')
            fiber_mesh =  bundles['aims_Tmtktri']
            aims.SurfaceManip.meshTransform(fiber_mesh.get(), 
                                            bundles_to_ROIs_motion)
            aims.GraphManip.storeAims(roi_graph, edge.get(), 
                                      'roi_mesh_junction', fiber_mesh)

