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


#----------------------------Imports-------------------------------------------


from __future__ import print_function

# soma module
from soma import aims
from soma.aims.meshSplit import meshSplit2

# anatomist
import anatomist.cpp as anatomist
import anatomist.direct.api as ana

# constel module
from constel.lib.utils.graphtools import merge_bundlegraph_and_roigraph
import constel


#----------------------------Functions-----------------------------------------


class BundlesSplitByCorticalROIsModule(anatomist.Module):
    def name(self):
        return "Bundles split by cortical regions module"

    def description(self):
        return ("Bundles split by cortical ROIs using a cortical mesh,"
                " labels texture, and bundles graph")


class FusionBundlesSplitByCorticalROIsMethod(anatomist.FusionMethod):

    def __init__(self):
        anatomist.FusionMethod.__init__(self)

    def canFusion(self, objects):
        # number of arguments
        if len(objects) not in (3, 4):
            return False

        mesh = None
        tex = None
        bundles_graph = None

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
                return False

        # fusion check
        if mesh and bundles_graph and tex:
            return True
        else:
            return False

    def fusion(self, objects):
        # type-checking
        ima = None
        for obj in objects:
            if isinstance(obj, anatomist.ASurface_3):
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
        if ima:
            aims_ima = anatomist.AObjectConverter.aims(ima)
        aims_mesh = anatomist.AObjectConverter.aims(mesh)

        # load an existing transformation
        dw_to_t1 = a.getTransformation(ref_bundles, ref_anat)
        if dw_to_t1:
            motion = dw_to_t1.motion()
        else:
            raise TypeError("Non valid transformation matrix: dw to t1")

        # 3D resolution
        voxel_size = None  # no voxels by default
        if ima:
            voxel_size = aims_ima.header()["voxel_size"]

        # handle texture timestep
        tex_hdr = aims_tex.header()
        if "time_step" in tex_hdr:
            time_step = tex_hdr["time_step"]
        else:
            time_step = 0

        # number of vertices
        n = aims_tex[time_step].nItem()

        if n != 0:
            filter_proportion = 0
            if "fibers_proportion_filter" in bundles_graph.graph():
                filter_proportion \
                    = bundles_graph.graph()["fibers_proportion_filter"]
            aims_bundles_graph = constel.texMeshAndBundles_to_BundlesGraph(
                aims_mesh, aims_tex, "Name1_Name2", bundles_graph.fileName(),
                motion, '', filter_proportion, time_step)

            # creating a roi graph: voxel size and mesh corresponding
            aims_roi_graph = aims.Graph("RoiArg")
            if voxel_size is not None:
                aims_roi_graph["voxel_size"] = voxel_size
            else:
                aims_roi_graph["voxel_size"] = [1., 1., 1.]
            aims_roi_graph["global_mesh"] = mesh

            # roi_mesh and aims_roi (bucket) are added in aims_roi_graph
            meshSplit2(aims_mesh,
                       aims_tex,
                       aims_roi_graph,
                       voxel_size=voxel_size,
                       tex_time_step=time_step)

            # merge Aims ROIs graph with Aims bundles graph (based on the same
            # transformation matrix)
            # Note: the "background" node is "other" in aims_bundles_graph
            #   and "0" in aims_roi_graph
            merge_bundlegraph_and_roigraph(aims_roi_graph,
                                           aims_bundles_graph,
                                           motion)
                                           #nodes_names_mapping={"others": 0})

            # converts an Aims object to an Anatomist object
            roi_graph = anatomist.AObjectConverter.anatomist(aims_roi_graph)

            # color for roi_mesh and roi_mesh_junction (label and filename)
            aims.GraphManip.setAttributeColor(aims_roi_graph,
                                              "roi_mesh",
                                              [128, 0, 255, 20])
            aims.GraphManip.setAttributeColor(aims_roi_graph,
                                              "roi_mesh_junction",
                                              [10, 0, 255, 128])
            aims_roi_graph["aims_reader_loaded_objects"] = 3

            meshes = []
            a = ana.Anatomist()
            for vertex in aims_roi_graph.vertices():
                if "roi_mesh_ana" in vertex:
                    meshes.append(a.AObject(a, vertex["roi_mesh_ana"]))
            if meshes:
                a.execute("SetMaterial",
                          objects=meshes,
                          selectable_mode="always_selectable")

            # property
            roi_graph.setName("ROIs And Bundles Graph")
            roi_graph.setReferential(ref_anat)
            roi_graph.updateAfterAimsChange()
            roi_graph.notifyObservers()
            del aims_bundles_graph
            del aims_roi_graph

        else:
            print("error in TimeTexture: 0 Items")

        # remove python references to anatomist objects before closing
        del bundles_graph, tex, obj
        return roi_graph

    def ID(self):
        return "FusionBundlesSplitByCorticalROIsMethod"
#FusionBundlesSplitByCorticalRegions

f = anatomist.FusionFactory.factory()
m = FusionBundlesSplitByCorticalROIsMethod()
f.registerMethod(m)

pm = BundlesSplitByCorticalROIsModule()
