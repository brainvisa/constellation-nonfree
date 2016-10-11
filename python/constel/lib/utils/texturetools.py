#!/usr/bin/env python
###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################

from __future__ import print_function
from soma import aims
import numpy as np

# System import
import logging
import numpy
import re

# Define logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Create formatter (time and level of each message)
formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')

# Create handler (console)
steam_handler = logging.StreamHandler()
steam_handler.setLevel(logging.INFO)
# steam_handler.setFormatter(formatter)
logger.addHandler(steam_handler)


def identify_patch_number(fname):
    """Identify the number of the patch from a filename.

    Parameter
    ---------
    fname: str (mandatory)
        the name of profile containing the patch number

    Return
    ------
    patch_nb: int
        the patch number
    """
    patch = re.findall("G[0-9]+", fname)[0]
    patch_nb = int(re.findall("[0-9]+", patch)[0])
    return patch_nb


def management_internal_connections(patch, mask, profile, list_regions=[]):
    """Remove or keep the patch internal connections.

    Parameters
    ----------
    patch: str (mandatory)
        the number of ROI
    mask: numpy.ndarray (mandatory)
        a labeling of gyri cortical segmentation
    profile: numpy.ndarray (mandatory)
        a cortical connectivity profile
    list_regions: list of int (mandatory)
        list of labels (int values) corresponding to gyri

    Return
    ------
    aims_profile: numpy.ndarray
        a cortical connectivity profile
        with or without the path internal connections
    """

    name_labels = list(numpy.unique(mask))
    for element in list(list_regions):
        if element in name_labels:
            name_labels.remove(element)
    for name_label in name_labels:
        profile[mask == name_label] = 0

    return profile


def remove_labels(tex, labels):
    """Allows to remove labels in a texture from a given labels list.

    Parameters
    ----------
    tex (aims time texture S16):
        labeled texture between 1 and max.
    labels (list):
        labels to remove

    Return
    ------
    otex (aims time texture S16): labeled texture,
                                  without region of labels.
                                  renumbered between 1 and new max.
                                  (max - len(labels))
    """
    # create the aims file
    otex = aims.TimeTexture_S16()
    tex_ar = tex[0].arraydata()
    otex[0].assign(tex_ar)
    otex_ar = otex[0].arraydata()

    # replace the labels to be deleted by zero (background)
    for l in labels:
        otex_ar[otex_ar == l] = 0

    # keep only the labels different of zero
    otex_kept_labels = np.unique(otex_ar)
    otex_kept_labels_list = otex_kept_labels.tolist()
    if otex_kept_labels_list.count(0) != 0:
        otex_kept_labels_list.remove(0)

    # relabelize all the labels from 1
    for i in xrange(len(otex_kept_labels_list)):
        current_label = otex_kept_labels_list[i]
        print("current label:", current_label, " new:", str(i + 1))
        otex_ar[tex_ar == current_label] = i + 1
    return otex


def texture_time(k_max, clusters_id, vertices_patch, vertices_mesh, mode,
                 minid=None, maxid=None):
    """Gives a parcellation per clusters number in a given region.

    Parameters
    ----------
    k_max (int):
        number of clusters.
    clusters_id (list):
        list of arrays. Each array is used to store the cluster number
        to which each item was assigned by the clustering algorithm.
    vertices_mesh (long):
        number of vertices for a given mesh.
    vertices_patch (list):
        set of vertices for a given patch.
    mode (int):
        1 -- average.
        2 -- concatenate.
    Return
    ------
    tex (aims time texture S16):
        for each time step, results of the k-th clustering.
    """
    tex = aims.TimeTexture_S16()
    for k in range(k_max - 1):
        tex[k].resize(vertices_mesh, 0)
        if mode == 1:
            tex[k].arraydata()[vertices_patch] = \
                clusters_id[k].astype(np.int16)
        if mode == 2:
            print(vertices_patch)
            tex[k].arraydata()[vertices_patch] = \
                clusters_id[k][minid:maxid].astype(np.int16)

    return tex


def add_texture(tex, add_tex, add_value):
    """Allows to add texture for another.

    Parameters
    ----------
    tex (aims time texture S16):
        patch texture 1, fixed time-step.
    add_tex (aims time texture S16):
        patch texture 2, fixed time-step.
    add_value (int):
        value added to the original numbering.

    Return
    ------
    tex: add_tex + add_value, keep its value otherwise.
    """
    vertices_tex = int(tex[0].nItem())
    vertices_added = int(add_tex[0].nItem())
    if vertices_tex != vertices_added:
        raise ValueError("tex and tex_add have not the same size: "
                         + str(vertices_tex) + " and "
                         + str(vertices_added))
    for v in xrange(vertices_tex):
        v_tex_add = add_tex[0][v]
        v_tex = tex[0][v]
        if v_tex_add > 0 and v_tex == 0:
            tex[0][v] = v_tex_add + add_value
        elif v_tex_add > 0 and v_tex > 0:
            raise ValueError("tex and add_tex have a common non zero value")
    return tex


def normalize_profile(profile):
    """Normalization of connectivity profile by its total number of
    connections.

    Parameter
    ---------
    profile: numpy.ndarray (mandatory)
       a cortical connectivity profile

    Return
    ------
    profile: numpy.ndarray
       normalized profile
    """
    sum_values = profile.sum()
    if sum_values > 0:
        z = 1. / sum_values
        profile *= z
    logger.info("The profile is normalized.")

    return profile


def geodesic_gravity_center(mesh, parcellated_texture, region_name,
                            time_step=0):
    """Compute the gravity center of the cluster or given region.

    Parameters
    ----------
    mesh: (AimsTimeSurface_3_VOID)
    parcellated_texture: (TimeTexture_S16)
    region_name: (int)
    time_step: (int)

    Return
    ------
    gravity_center_index: (int)

    .. note::
        (Approximated)
        For circular geodesic region it works, otherwise the max
        of distance_texture is not single.
    """

    # get the number of vertices of the mesh
    number_of_vertices = mesh.vertex().size()

    # create a texture of type TimeTexture_S16 (soma.aims)
    # with an equal amount of vertices as the mesh
    texture = aims.TimeTexture_S16()
    texture[0].resize(number_of_vertices, 1)

    # get a list of vertices linked to the given region, not the entire mesh
    vertices_index = numpy.where(
        parcellated_texture[time_step].arraydata() == region_name)[0]

    # each elements of vertices_index is remplaced by 0,
    # the others are at 1 in texture
    for (index, vertex) in enumerate(vertices_index):
        texture[0][vertex] = 0

    # compute geodesic depth of a triangulation,
    # the texture define the object and the background,
    # the radius-*are morphological parameter*
    distance_texture = aims.meshdistance.MeshDistance(mesh, texture, False)

    # get the max to define the gravity center of the given region
    gravity_center_index = distance_texture[0].arraydata().argmax()

    return gravity_center_index


def concatenate_texture(cortical_parcellations, time_step):
    """Concatenate all the cortical parcellations on the same mesh.

    Parameters
    ----------
    cortical_parcellations: list of str (mandatory)
        List of the cortical parcellations of a subject.
    timestep: list of int (mandatory)
        List of all time step for the cortical parcellations.

    Returns
    -------
    final_rseg: str
        The complete cortical parcellation.
    """
    for idx, filename in enumerate(cortical_parcellations):
        roiseg = aims.read(filename)
        rseg = numpy.array(roiseg[time_step[idx]].arraydata())
        if idx == 0:
            tmp_rseg = rseg
            max_time_step = max(rseg)
        else:
            l = numpy.unique(rseg)
            if l[0] == 0:
                labels = l.nonzero()[0][::-1]
            else:
                labels = l[::-1]
            for label in labels:
                rseg[rseg == label] = label + max_time_step
            temp = tmp_rseg[:]
            tmp_rseg = [x + y for x, y in zip(temp, rseg)]
            max_time_step = max(rseg)
    final_rseg = aims.TimeTexture_S16()
    final_rseg[0].assign(tmp_rseg)
    return final_rseg


def create_relationship_region2neighbors(meshname, segname):
    """
    """
    # load the files (mesh and gyri segmentation)
    mesh = aims.read(meshname)
    tex = aims.read(segname)

    # create a numpy array from gyri segmentation (aims to numpy)
    texar = tex[0].arraydata()

    # load the polygons of the mesh (numpy array)
    triangles = numpy.asarray([tri.arraydata() for tri in mesh.polygon()])

    dict_neighboors = {}
    for search in xrange(len(texar)):
        neighboors = []

        # give the polygons having the nodes "search"
        idx = numpy.where(triangles == search)

        # give the index of the nodes
        for d in idx[0]:
            neighboors.extend(triangles[d])

        # list of the neighboors and give the labels
        neighboors = numpy.unique(numpy.asarray(neighboors))
        neighboors = neighboors[neighboors != search]
        labels = []
        for n in neighboors:
            labels.append(texar[n])
        labels = numpy.unique(labels)
        labels = labels[labels != texar[search]]

        # complete the dictionnaire with the new labels
        if texar[search] in dict_neighboors:
            labelsindico = dict_neighboors[texar[search]]
            for l in labels:
                if l not in labelsindico:
                    dict_neighboors[texar[search]].append(l)
        else:
            dict_neighboors[texar[search]] = labels.tolist()
    for i in dict_neighboors.values():
        i.sort()

    mintex = min(texar)
    maxtex = max(texar)

    # load the labels
    countries = range(mintex, maxtex + 2)

    for i in countries:
        if i not in dict_neighboors.keys():
            if i != (len(dict_neighboors.keys()) + 1):
                dict_neighboors[i] = []

    return dict_neighboors
