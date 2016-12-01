#!/usr/bin/env python

# python system modules
from __future__ import print_function
from tempfile import mkdtemp
import os.path as path
import numpy
import struct
import os

# scipy
import scipy
import scipy.cluster.vq
import scipy.spatial.distance
dst = scipy.spatial.distance.euclidean
from scipy.cluster.hierarchy import fcluster

# fastcluster
import fastcluster


def nearest_neighbour_profiles(profile, atlas, atlas_classes, indices):
    """ Compute the nearest neighbour profiles.
    
    Parameters
    ----------
    profile: numpy.array (mandatory)
        Individual reduced matrix, shape (region vertices, basins).
    atlas: numpy.array (mandatory)
        Group reduced matrix, shape (region vertices, basins).
    atlas_classes: nupy.array (mandatory)
        Clustering of the cortical region.
    indices: numpy.array (mandatory)
        Vertex indices belonging to the cortical region.

    Returns
    -------
    res_classes: numpy.array
        Individual clustering of the cortical region defined from atlas.
    """
    res_classes = numpy.zeros(len(atlas_classes))
    for l, line in enumerate(profile):
        distance = numpy.sum((line - atlas)**2, axis=1)
        atlas_line = numpy.argmin(distance)
        i = indices[atlas_line]
        res_classes[indices[l]] = atlas_classes[i]
    return res_classes


def get_centers(clusters_id, kmax_clusters):
    """Get all clusters centers.

    Parameters
    ----------
    clusters_id: (array)
        From 1 to kmax_clusters, a label is attributed for each vertex,
        defining clusters
    kmax_clusters: (integer)
        the number of clusters to be found

    Return
    ------
    centers: (array)
        provides an array of floats corresponding to a center for each cluster
    """
    i = 0
    centers = numpy.array([])
    while (centers.size < kmax_clusters):
        print(i)
        c = clusters_id[i]
        if ((numpy.where(centers == c)[0]).size == 0):
            centers = numpy.append(centers, numpy.array([int(c)]))
        i += 1
    return centers


def get_centroid(distance_matrix):
    """Get the centroid of all data.

    It is the sample that minimizes the distance to all other samples
    therefore it is the sample with the smallest sum of values on its column.

    Parameter
    ---------
    distance_matrix: ()
        The distance matrix is symmetric and has zeros on the diagonal

    Return
    ------
    centroid: (integer)
        value of the centroid
    """
    # sum columns of the matriwcDistx
    s = distance_matrix.sum(axis=0)
    centroid = numpy.where(s == s.min())[0][0]
    return centroid


def get_intradistance(distance_matrix, clusters_id, kmax_clusters, factor):
    """Get within cluster sum of distances.

    Parameters
    ----------
    distance_matrix: ()
        The distance matrix is symmetric and has zeros on the diagonal
    clusters_id: (array)
        From 1 to kmax_clusters, a label is attributed for each vertex,
        defining clusters
    kmax_clusters: (integer)
        the number of clusters to be found
    factor: (integer)
        in generally, factor = 2 for the square of the distances

    Return
    ------
    addW: (interger)
        within cluster sum of distances
    """
    centers = get_centers(clusters_id, kmax_clusters)
    W = numpy.zeros(kmax_clusters)
    addW = 0.0
    for l in range(kmax_clusters):
        sum_distance = 0.0
        centerid = numpy.where(clusters_id == centers[l])[0]
        for i in centerid:
            for j in centerid:
                if (i < j):
                    sum_distance += (distance_matrix[j][i] ** factor)
                elif (j < i):
                    sum_distance += (distance_matrix[i][j] ** factor)
        W[l] = sum_distance / float(2.0 * centerid.size)
        addW += W[l]
    return addW


def distToCenters(distance_matrix, centers, clusters_id, kmax_clusters):
    """Get the average distance to the centroid, for each cluster.

    Parameters
    ----------
    distance_matrix: ()
        The distance matrix is symmetric and has zeros on the diagonal
    centers: (arrauy)
        provides an array of floats corresponding to a center for each cluster
    clusters_id: (array)
        From 1 to kmax_clusters, a label is attributed for each vertex,
        defining clusters
    kmax_clusters: (integer)
        the number of clusters to be found

    Return
    ------
    dist_to_center: ()
        ...
    """
    i = 0
    dist_to_center = numpy.zeros(kmax_clusters)
    for i in range(kmax_clusters):
        c = centers[i]
        group = numpy.where(clusters_id == c)[0]
        dgroup = distance_matrix[c, group]
        dist_to_center[i] = dgroup.mean()
    return dist_to_center


def sse(tab):
    """Sum square error function"""
    v = 0.5*numpy.sum(numpy.square(tab))
    return v


def entropy(labels):
    '''Computes entropy of label distribution'''
    n_labels = len(labels)
    if n_labels <= 1:
        return 0
    counts = numpy.bincount(labels)
    counts = numpy.array([x for x in counts if x != 0]).astype(float)
    counts_sum = numpy.sum(counts)
    probs = counts / n_labels
    n_classes = numpy.count_nonzero(probs)
    if n_classes <= 1:
        return 0
    entropy = -numpy.sum(probs * (numpy.log(counts) - numpy.log(counts_sum)))
    return entropy


def gap(data, refs=None, nrefs=20, ks=range(1, 11)):
    """Compute the Gap statistic for an nxm dataset in data.

    Either give a precomputed set of reference distributions in refs as an
    (n,m,k) scipy array, or state the number k of reference distributions
    in nrefs for automatic generation with a uniformed distribution within
    the bounding box of data.

    Give the list of k-values for which you want to compute the statistic
    in ks.
    """
    shape = data.shape
    if refs is None:
        tops = data.max(axis=0)
        bots = data.min(axis=0)
        dists = scipy.matrix(scipy.diag(tops-bots))

        rands = scipy.random.random_sample(size=(shape[0], shape[1], nrefs))
        for i in range(nrefs):
            rands[:, :, i] = rands[:, :, i] * dists + bots
    else:
        rands = refs

    gaps = scipy.zeros((len(ks),))
    for (i, k) in enumerate(ks):
        print("ks", ks, "k", k, "i", i)
        # clusterid, error, nfound = pc.kmedoids(data, k, 100)
        # labels, error, nfound = pc.kcluster(clusterid, k)
        # kml = labels
        # centers = get_centers(clusterid, k)
        # kmc = centers
        (kmc, kml) = scipy.cluster.vq.kmeans2(data, k)
        print("labels -->", kml, "centers -->", kmc)
        disp = sum([dst(data[m, :], kmc[kml[m], :]) for m in range(shape[0])])

        refdisps = scipy.zeros((rands.shape[2],))
        for j in range(rands.shape[2]):
            (kmc, kml) = scipy.cluster.vq.kmeans2(rands[:, :, j], k)
            refdisps[j] = sum(
                [dst(rands[m, :, j], kmc[kml[m], :]) for m in range(shape[0])])
        gaps[i] = scipy.log(scipy.mean(refdisps)) - scipy.log(disp)

    return gaps


def ward_method(dmat_file, n, output_dir, n_clusters):
    """Performs Wards linkage on the distance matrix.

    Ward's method is a criterion applied in hierarchical cluster analysis.
    Ward's minimum variance criterion minimizes the total within-cluster
    variance.

    For disjoint clusters C_i, C_j, and C_k with sizes n_i, n_j, and n_k.

                             (n_i + n_k)
    D(C_i U C_j, C_k) =  ----------------- D(C_i, C_k) +
                         (n_i + n_j + n_k)
                            (n_i + n_k)
                         ----------------- D(C_j, C_k) -
                         (n_i + n_j + n_k)
                                n_k
                         ----------------- D(C_i, C_j)
                         (n_i + n_j + n_k)
    Parameters
    ----------
        dmat_file (array[..., ]): binary file of the distance matrix
        n (int): shape[0] of the reduced matrix
        output_dir (str): directory to save the dendrogram
        n_clusters (int): number of clusters (kmax) for the clustering

    Returns
    -------
        clusterid (array[..., ]):
                  From 1 to kmax, a label is attributed for each vertex,
                  defining clusters.

    **References**
        .. [1] Ward, J. H.
               Hierarchical grouping to optimize an objective function.
               Journal of the American Statistical Association,
               58:238-244, 1963.
        .. [2] Wishart, D.
               An algorithm for hierarchical classification.
               Biometrics 25:165-170, 1969.
        .. [3] Cormack, R. M.
               A Review of Classification.
               Journal of the Royal Statistical Society,
               Series A, 134(3), 321-367, 1971.
    """

    # distance matrix file loading
    f = open(dmat_file, "r")

    # initialization of the distance matrix as vector
    # memory map
#    a = numpy.memmap('test.mymemmap', dtype=numpy.single, mode='w+', shape=(n*(n-1)/2,))
#    del a
#    dist_mat = numpy.memmap('test.mymemmap', dtype=numpy.single, mode='r+', shape=(n*(n-1)/2,))
    dist_mat = numpy.zeros((n*(n-1)/2,))

    # number of iteration to generate the distance matrix
    n_iter = 100000  # not perfect
    nb_iteration = dist_mat.shape[0] / n_iter
    nb_add = dist_mat.shape[0] % n_iter
    list_iteration = [nb_iteration] * n_iter

    # generate the distance matrix from dmat_file to dist_mat (array)
    idx = 0
    for j in xrange(n_iter):
        if j < nb_add:
            list_iteration[j] += 1
        x = f.read(list_iteration[j] * 8)
        dist_mat.ravel()[
            idx:idx + list_iteration[j]] = struct.unpack(
                'd' * list_iteration[j], x)
        idx += list_iteration[j]
    f.close()

    # save the dist_mat
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_distmat = output_dir + '/dist_mat.npy'
    numpy.save(output_distmat, dist_mat)

    # compute linkage ward
    Z = fastcluster.linkage(dist_mat, method='ward', preserve_input=False)

    # save the dendrogram
    output_Z = output_dir + '/dendrogram.npy'
    numpy.save(output_Z, Z)

    # labelization between 1 to nb
    clusterid = []
    for nb in range(1, n_clusters + 1):
        print("Trying {nb} cluster(s)".format(nb=nb))
        clusters = fcluster(Z, criterion='maxclust', t=nb)
        clusterid.append(clusters)
    return clusterid
    
    # save the clusterid
    output_clusterid = output_dir + '/clusterid.npy'
    numpy.save(output_clusterid, clusterid)


def calinski_harabasz_index(distance_matrix, cluster_id, kmax_clusters):
    """Evaluate the optimal number of clusters.

    Calinski and Harabasz propose to choose the value of K which maximizes an
    index based on the quotient:

                             [trace B / (K - 1)]
                     CH(K) = -------------------   for K E N
                             [trace W / (n - K)]

    where B denotes the error sum of squares between different clusters
    (inter-cluster) and W the squared differences of all objects in a cluster
    from their respective cluster center (intra-cluster).

    Calculated for each possible cluster solution the maximal achieved index
    value indicates the best clustering of the data. The important
    characteristic of the index is the fact that on the one hand trace W will
    start at a comparably large value. With increasing number of clusters K,
    approaching the optimal clustering solution in K* groups, the value
    should significantly decrease due to an increasing compactness of each
    cluster. As soon as the optimal solution is exceeded an increase in
    compactness and thereby a decrease in value might still occur; this
    decrease, however, should be notably smaller. On the other hand, trace T
    should behave in the opposite direction, getting higher as the number of
    clusters K increases, but should also reveal a kind of softening in its
    rise if K gets larger than K*.

    Parameters
    ----------
    distance_matrix:
        The distance matrix is symmetric and has zeros on the diagonal
    cluster_id:
        This array stores the cluster number to which each item was assigned
        by the clustering algorithm
    kmax_clusters: (int)
        the number of clusters to be found

    Return
    ------
    ch_index:
        index value indicates the best clustering of the data.
        maximum value of the index
    """
    # initialization of index value and trace B, W
    ch_index = 0
    B = 0
    W = 0

    centers = get_centers(cluster_id, kmax_clusters)
    centroids = get_centroid(distance_matrix)

    n_sample = distance_matrix.shape[0]

    for i in range(kmax_clusters):
        # trace B
        B += (distance_matrix[centers[i], centroids] * numpy.where(
            cluster_id == centers[i])[0].size)
        for j in range(n_sample):
            if cluster_id[i] == centers[j]:
                # trace W
                W += distance_matrix[j, centers[i]]
    # index values
    ch_index = B * (n_sample - kmax_clusters) / (
        W * 10.0 * (kmax_clusters - 1))
    return ch_index


def davies_bouldin_index(distance_matrix, cluster_id, kmax_clusters):
    """Evaluate the optimal number of clusters.

    Instead of simply proposing a cluster index, Davies and Bouldin formulated
    a general framework for the evaluation of the outcomes of cluster
    algorithms:

                             1   K
                    DB(K) =  -  sum(R_k)  for K E N
                             K  k=1

    For each cluster C_k an utmost similar cluster- regarding their
    intra-cluster error sum of squares is searched, leading to R_k. The index
    the defines the average over these values. In contrast to the
    aforementioned cluster indexes, here, the minimal observed index indicates
    the best cluster solution.

    Parameters
    ----------
    distance_matrix:
        The distance matrix is symmetric and has zeros on the diagonal
    cluster_id: (array)
        This array stores the cluster number to which each item was assigned
        by the clustering algorithm
    kmax_clusters: (int)
        the number of clusters to be found

    Return
    ------
    db_index:
        the minimal observed index indicates the best cluster solution
        minimum value of the index
    """
    db_index = 0

    centers = get_centers(cluster_id, kmax_clusters)
    distC = distToCenters(distance_matrix, centers, cluster_id, kmax_clusters)

    for i in range(kmax_clusters):
        Di = numpy.zeros(kmax_clusters)
        for j in range(kmax_clusters):
            if i != j:
                Di[j] = (distC[i] + distC[j]) / distance_matrix[
                    centers[i], centers[j]]
        db_index += Di.max()
    db_index = db_index / float(kmax_clusters)

    return db_index


def dunn_index(distance_matrix, clusters_id):
    """Evaluate the optimal number of clusters.

    Parameters
    ----------
    distance_matrix: ()
        The distance matrix is symmetric and has zeros on the diagonal
    clusters_id: (array)
        This array stores the cluster number to which each item was assigned
        by the clustering algorithm

    Return
    ------
    dunn_index: ()
        the minimal observed index indicates the best cluster solution
        maximum value of the index
    """
    distance_intra = numpy.array([])
    distance_inter = numpy.array([])

    n_sample = distance_matrix.shape[0]

    for i in range(n_sample):
        for j in range(i + 1, n_sample):
            if (clusters_id[i] != clusters_id[j]):
                distance_inter = numpy.append(
                    distance_inter, distance_matrix[i, j])
            else:
                distance_intra = numpy.append(
                    distance_intra, distance_matrix[i, j])
    dunn_index = distance_inter.min() / distance_intra.max()
    return dunn_index

def krzanowski_lai_index(distance_matrix, clusters_id, kmax_clusters):
    """Evaluate the optimal number of clusters.

    Krzanowski and Lai developed a cluster index that, similar to index of
    Calinski and Harabasz, is based on the squared differences o f all objects
    in a cluster from their respective cluster center - trace W.
    The authors define DIFF(K) as the difference between a clustering of the
    data in K and a clustering in K - 1 clusters. Let J be the number of
    variables that has been measured on each x_i E X and trace W_k the sum of
    squares function that corresponds to the clustering in K clusters.

           DIFF(K) = (K-1)^(2/J) . traceW_(K-1) - K^(2/J) . traceW_K

    The authors claim that if there exists an optimal clustering solution in
    K* groups, the value of DIFF(K) should be comparably large and positive.
    In contrast, all values of DIFF(K) for K > K* will have rather small
    values (maybe even negative), while values for K < K* will ba rather large
    and positive.

                                   DIFF(K)
                       KL(K) = | ----------- |
                                 DIFF(K + 1)

    Parameters
    ----------
    distance_matrix: (asarray)
        The distance matrix is symmetric and has zeros on the diagonal
    kmax_clusters: (int)
        th optimal number of clusters to be found
    clusters_id: (array)
        This array stores the cluster number to which each item was assigned
        by the clustering algorithm
    
    Return
    ------
    kl_index:
        the optimal cluster solution is indicated by the highest value of
        KL(K)
        maximum value of the index
    """
    kl_index = 0
    W = 0    

    centers = get_centers(clusters_id, kmax_clusters)
    centroids = get_centroid(distance_matrix)

    
