# -*- coding: utf-8 -*-
import numpy as np
from soma import aims
import anatomist.cpp as anatomist
from soma.aims.meshSplit import meshSplit2
from constel.lib.graphtools import mergeBundlesGraphAndROIsGraph
import constel


class MeshFusionMeshRoiGraphModule(anatomist.Module):
    def name(self):
        return 'Mesh RoiGraph Fusion Module'

    def description(self):
        return 'Merge a Mesh and a ROI Graph and create the associated Graph'


class FusionTexMeshImaAndBundlesToROIsAndBundlesGraphMethod(
        anatomist.FusionMethod):
    def __init__(self):
        anatomist.FusionMethod.__init__(self)

    def canFusion(self, objects):
        # number of arguments 
        if len(objects) != 4:
            raise TypeError('Fusion Method takes at least 4 arguments')

        # type-checking
        for obj in objects:
            if isinstance(obj, anatomist.ASurface_3):
                mesh = obj
            elif isinstance(obj, anatomist.AGraph):
                bundles_graph = obj
            elif obj.type() == anatomist.AObject.TEXTURE:
                tex = obj
            elif obj.type() == anatomist.AObject.VOLUME:
                ima = obj
            else:
                raise TypeError('Objects must be an Anatomist instance')

        # fusion check
        if (mesh and bundles_graph) and tex:
            return True
        else:
            raise TypeError('Could not fusion objects')

    def fusion(self, objects):
        # type-checking
        for obj in objects:
            if isinstance(obj, anatomist.ASurface_3 ):
                mesh = obj
            elif isinstance(obj, anatomist.AGraph):
                bundles_graph = obj
            elif obj.type() == anatomist.AObject.TEXTURE:
                tex = obj
            elif obj.type() == anatomist.AObject.VOLUME:
                ima = obj

        # instance of Anatomist
        a = anatomist.Anatomist()

        # get referential
        ref_bundles = bundles_graph.getReferential()
        ref_anat = mesh.getReferential()

        # converts an Anatomist object to an AIMS object
        aims_tex = anatomist.AObjectConverter.aims(tex)
        aims_ima = anatomist.AObjectConverter.aims(ima)
        aims_mesh = anatomist.AObjectConverter.aims(mesh)

        # load an existing transformation
        dw_to_t1 = a.getTransformation(ref_bundles, ref_anat)
        if dw_to_t1:
            motion = dw_to_t1.motion()
        else:
            raise TypeError('Non valid transformation matrix: dw to t1')

        # 3D resolution
        voxel_size = aims_ima.header()['voxel_size']

        # number of vertices
        n = int(aims_tex.nItem())

        if n != 0:
            filter_proportion = 0
            if bundles_graph.graph().has_key('fibers_proportion_filter'):
                filter_proportion \
                    = bundles_graph.graph()['fibers_proportion_filter']
            aims_bundles_graph = constel.texMeshAndBundles_to_BundlesGraph(
                aims_mesh, aims_tex, "Name1_Name2", bundles_graph.fileName(),
                motion, '', filter_proportion)

            # creating a roi graph: voxel size and mesh corresponding
            aims_roi_graph = aims.Graph('RoiArg')
            aims_roi_graph['voxel_size'] = voxel_size
            aims_roi_graph['global_mesh'] = mesh

            # roi_mesh and aims_roi (bucket) are added in aims_roi_graph
            meshSplit2(aims_mesh, aims_tex, aims_roi_graph, voxel_size)

            # merge Aims ROIs graph with Aims bundles graph (based on the same transformation matrix)
            mergeBundlesGraphAndROIsGraph(aims_roi_graph, aims_bundles_graph, 
                motion)

            # converts an Aims object to an Anatomist object
            roi_graph = anatomist.AObjectConverter.anatomist(aims_roi_graph)

            # color for roi_mesh and roi_mesh_junction (label and filename)
            aims.GraphManip.setAttributeColor(aims_roi_graph, 'roi_mesh',
                [128, 0, 255, 20])
            aims.GraphManip.setAttributeColor(aims_roi_graph, 
                'roi_mesh_junction', [10, 0, 255, 128])
            aims_roi_graph["aims_reader_loaded_objects"] = 3

            # property
            roi_graph.setName("ROIs And Bundles Graph")
            roi_graph.setReferential(ref_anat)
            roi_graph.updateAfterAimsChange()
            roi_graph.notifyObservers()
            del aims_bundles_graph
            del aims_roi_graph

        else:
            print "error in TimeTexture: 0 Items"


        # remove python references to anatomist objects before closing
        del bundles_graph, tex, obj
        return roi_graph

    def ID(self):
        return "FusionTexMeshImaAndBundlesToROIsAndBundlesGraphMethod"


f = anatomist.FusionFactory.factory()
m = FusionTexMeshImaAndBundlesToROIsAndBundlesGraphMethod()
f.registerMethod(m)

pm = MeshFusionMeshRoiGraphModule()
