#!/usr/bin/env python
import soma.aims.texturetools as MTT
from soma import aims
import numpy as np
import numpy


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
    otex = aims.TimeTexture_S16()
    tex_ar = tex[0].arraydata()
    otex[0].assign(tex_ar)
    otex_ar = otex[0].arraydata()
    for l in labels:
        otex_ar[otex_ar == l] = 0
    otex_kept_labels = np.unique(otex_ar)
    otex_kept_labels_list = otex_kept_labels.tolist()
    if otex_kept_labels_list.count(0) != 0:
        otex_kept_labels_list.remove(0)
    for i in xrange(len(otex_kept_labels_list)):
        current_label = otex_kept_labels_list[i]
        print "current label:", current_label, " new:", str(i + 1)
        otex_ar[tex_ar == current_label] = i + 1
    return otex


def change_wrong_labels(cc_label, label, gyri_tex, mesh_neighbors_vector,
                        cc_tex_label):
    """After a study of its neighbors, wrong label is replaced by the
    correct number.

    Parameters
    ----------
    cc_label: label of connected component in cc_tex_label
    label: label of associated vertices in gyri texture
    gyri_tex (aims time texture S16): gyri texture
    mesh_neighbors_vector : aims.SurfaceManip.surfaceNeighbours(mesh)
    cc_tex_label : texture representing connected components of label

    Returns
    -------
    gyri_tex (aims time texture S16): new gyri_tex texture,
                                      without isolated vertex.
    winner_label: the correct number.
    """
    indexes = np.where(cc_tex_label == cc_label)[0]
    neighbor_labels = []
    print 'Nb of wrong indexes: ', indexes.size
    for i in indexes:
        for n in mesh_neighbors_vector[i]:
            n_label = gyri_tex[0][n]
            if n_label != label:
                neighbor_labels.append(n_label)
    v_labels = np.unique(neighbor_labels)
    max_count = 0
    winnner_label = -1
    for l in v_labels:
        nb_v_labels = neighbor_labels.count(l)
        if nb_v_labels > max_count:
            print 'Number of neighbor labels: ', nb_v_labels, 'for', l
            winnner_label = l
            max_count = nb_v_labels
    for i in indexes:
        gyri_tex[0][i] = winnner_label
    return gyri_tex, winnner_label


def find_wrong_labels(mesh, gyriTex):
    """

    Parameters
    ----------
    mesh:
    gyriTex: gyri texture

    Returns
    -------
    wrong_labels: list of wrong labels
        [cctex: connectedComponentTex: aimsTimeTexture_S16,
        time step = LabelsNb, for each time step (label in the tex),
        texture of the connected components corresponding to this label
        (background = -1, and connected components = values between 1 and ccNb)
        areas_measures = python dictionary,
        areas_measures[label] = [16.5, 6.0] (numpy array)
        if label (in tex) has two connected Components 1 and 2 with
        area = 16.5 and 6.0 respectively, areas are in square mm]
    """

    meshNeighborsVector = aims.SurfaceManip.surfaceNeighbours(mesh)
    cctex, areas_measures = MTT.connectedComponents(mesh, gyriTex, areas_mode=1)
    wrong_labels = []
    for label in areas_measures.keys():
        if areas_measures[label].size != 1:
            wrong_labels.append(label)
    return wrong_labels


def clean_gyri_texture(mesh, gyri_tex):
    """Cleaning a gyri texture by using connected components.

    Parameters
    ----------
    mesh (aims time surface):
        white mesh associated to gyri_tex
    gyri_tex (aims time texture S16):
        gyri texture as full FreeSurfer parcellation.
    Return
    ------
    gyri_tex (aims time texture S16):
        new gyri texture, without isolated vertex.
    """
    mesh_neighbors_vector = aims.SurfaceManip.surfaceNeighbours(mesh)
    cc_tex, areas_measures = MTT.connectedComponents(
        mesh, gyri_tex, areas_mode=1)
    wrong_labels = []
    for label in areas_measures.keys():
        if areas_measures[label].size != 1:
            wrong_labels.append(label)
    for label in wrong_labels:
        cc_tex_label = cc_tex[label - 1].arraydata()
        areas_measures_cc = areas_measures[label]
        cc_nb = areas_measures_cc.size
        for l in range(1, cc_nb):
            gyri_tex, win = change_wrong_labels(
                l + 1, label, gyri_tex, mesh_neighbors_vector, cc_tex_label)
    return gyri_tex


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


def normalized_texture(texture):
    """
    """
    tex = texture[0].arraydata()
    dividende_coef = tex.max()
    if dividende_coef > 0:
        z = 1./dividende_coef
        for i in xrange(texture[0].nItem()):
            value = texture[0][i]
            tex[0][i] = z*value
    return texture


def geodesic_gravity_center(mesh, parcellated_texture, region_name,
                            time_step=0):
    """Compute the gravity center of the cluster or given region.

    Parameters
    ----------
    mesh: (AimsTimeSurface_3)
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