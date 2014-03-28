#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from math import log
from scipy.misc import comb
import constel.lib.connmatrix.connmatrixtools as clccm

def intersection(list1, list2):
    """Allows to calculate the cardinal of the intersection of two regions.
    
    Args: 
        sets of indexes corresponding to 2 input regions (lists):
            list1
            list2
    Return:
        cardinal of the intersection.

    """
    intersection = 0
    for i in list2:
        intersection += list1.count(i)
    return intersection
  
def union(list1, list2):
    """Allows to calculate the cardianl of the union of two regions.
    
    Args:
        sets of indexes corresponding to 2 input regions (lists):
            list1
            list2
    Return:
        cardinal of the union.

    """
    intersection = intersection(list1, list2)
    union = len(list1) + len(list2) - intersection
    return union

def mutual_information(list1, list2):
    """The mutual information is a measure of the variables mutual dependence.
    
    Args:
        list1:
        list2:
    Return:
        mi (float): mutual information
    
    """
    contingency_matrix = clccm.contingency_matrix(list1, list2)
    contingency = np.array(contingency_matrix, dtype = 'float') 
    contingency_sum = np.sum(contingency)
    pi = np.sum(contingency, axis = 1)
    pj = np.sum(contingency, axis = 0)
    outer = np.outer(pi, pj)
    nnz = contingency != 0.0
    contingency_nm = contingency[nnz]
    log_contingency_nm = np.log(contingency_nm)
    contingency_nm /= contingency_sum
    log_outer = -np.log(outer[nnz]) + log(pi.sum()) + log(pj.sum())
    mi = (contingency_nm * (log_contingency_nm - log(contingency_sum))
           + contingency_nm * log_outer)
    mi = mi.sum()
    return mi
  
def rand_index(list1, list2):
    """The Rand index is a measure of the similarity between two data clusterings.
    
    Args:
        list1:
        list2:
    Return:
        rand_id (float): rand index
    """
    contingency_matrix = clccm.contingency_matrix(list1, list2)
    Nsample = len(list1)
    sum_c1 = sum(comb( n_c, 2, exact = 1) for n_c in contingency_matrix.sum(axis = 1))
    sum_c2 = sum(comb( n_k, 2, exact = 1) for n_k in contingency_matrix.sum(axis = 0))
    sum_ = sum(comb(nij, 2, exact = 1) for nij in contingency_matrix.flatten())
    prod_ = (sum_c1 * sum_c2) / float( comb( Nsample, 2 ) )
    mean_ = (sum_c2 + sum_c1) / 2.
    rand_id = ((sum_ - prod_) / (mean_ - prod_))
    return rand_id

def bohlandIndex(list1, list2):
    """Probability that an element is in set1 given that it is on in set2
    assymetric. Bohland, "The brain Atlas Concordance Problem", PlosOne, 2009
    
    Args:
        sets of indexes corresponding to to input regions, type: lists
            list1
            list2
    Return:
        bohland_id (float): bohland index between 0 and 1.
    
    """
    card_set2 = len(list2)
    if card_set2 <= 0:
        raise exceptions.ValueError("set2 is empty")
    intersection = 0
    for i in list2:
        intersection += list1.count(i)
    bohland_id = float(intersection) / card_set2
    return bohland_id
  
def bohlandSymmetricIndex(list1, list2):
    """Symmetrized bohland index.
    
    Args:
        sets of indexes corresponding to to input regions, type: lists
            list1
            list2
    Return:
        bohland_sym_id (float): Symmetrized bohland index.
    """
    p12 = bohlandIndex(list1, list2)
    p21 = bohlandIndex(list2, list1)
    bohland_sym_id = sqrt(p12 * p21)
    return bohland_sym_id
    
def jacard_index(list1, list2):
    """ The Jaccard coefficient measures similarity between finite sample sets,
    and is defined as the size of the intersection divided by the size of 
    the union of the sample sets.
    
    Args:
        sets of indexes corresponding to two input regions (lists):
            list1
            list2
    Return:
        Jaccard index (float)
        
    """
    intersection = intersection(list1, list2)
    union = union(list1, list2)
    if union <= 0:
        raise exceptions.ValueError("both sets are empty")
    jaccard_id = float(intersection) / union
    return jaccard_id
  
def jaccard_distance(list1, list2):
    """The Jaccard distance measures dissimilarity between sample sets, 
    by dividing the difference of the sizes of the union and the intersection 
    of two sets by the size of the union
    
    Args:
        sets of indexes corresponding to two input regions (lists):
            list1
            list2
    Return:
        Jaccard distance (float)
    """
    jaccard_id = jacard_index(list1, list2)
    jaccard_dist = 1 - jaccard_id
    return jaccard_dist

def cramer_v(list1, list2):
    """Cramer's V is a measure of association between two nominal variables 
    and it is calculated based on chi-square statistic. Use the Cramer’s V 
    statistic to assess the relative strength of the derived association.
    Gives values within the interval [0, 1], 1 indicating a perfect match.
    
    Parameters:
    ----------
        - list1: clustering for a group of subject
        - list2: clustering for an other group of subject
    
    Returns:
    ------
        - Cramer's value
    """
    
    # background is removed
    l1 = [i for i in list1 if i != 0]
    l2 = [i for i in list2 if i != 0]
    
    # contingency table
    matrix = clccm.contingency_matrix(l1, l2) 
    
    # depend on sample size
    col_sum = matrix.sum(0)
    row_sum = matrix.sum(1)
    n = matrix.sum()
    
    # matrix dimension
    x, y = matrix.shape
    
    # k is the number of rows or the number of columns, whichever is less
    k = min(x, y)
    
    # degrees of freedom
    df = (x - 1) * (y - 1)
    
    # expected frequency of occurrence in each cell:
    # the probability that a cluster number represents a given class is 
    # given by the cluster’s proportion of the row total.
    ef = (row_sum * col_sum) / float(n)
    
    # chi square stats:
    # determine whether an association exists
    # do not measure the strength of the association
    chi2 = (np.array((matrix - ef)) ** 2) / ef
    chi2 = chi2.sum()
    
    # Cramer's V:
    # Measuring strength of an association, is independent of the order of 
    # cluster labelling in the cross-matching of one cluster allocation 
    # against another
    # Cramer’s V ranges from -1 to 1 for 2X2 tables
    es = np.sqrt(chi2 / (n * (k - 1)))
    
    return es    
    
def dunn_index():
  pass

def db_index(distance_matrix, labels):
  """
  Davies-Bouldin clustering evaluation index.
  matrix: n*dim
  labels: cluster numbers corresponding to matrix sample n*1
  """
  pass

def calinsky_harabasz_index():
  pass

def gap_index(matrix, b):
  """
  Gap clustering evaluation index.
  data: samples along the rows
  clusters: vector containing the number of clusters
  """
  matrix = matrix.transpose()
  print 'Opening and reading: ', matrix
  Nsample = matrix.shape[0]
  x, y = matrix[:,0], matrix[:,1]
  xmin = x.min()
  xmax = x.max()
  ymin = y.min()
  ymax = y.max()
  for i in xrange(b):
    Xb = np.vstack([np.random.uniform(xmin, xmax, Nsample), 
                    np.random.uniform(ymin, ymax, Nsample)]).T
    clusterid, error, nfound = pc.kmedoids(matrix, K, 100)
    
  pass