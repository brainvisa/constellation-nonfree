#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from math import log
from scipy.misc import comb
import constel.lib.connmatrix.connmatrixtools as clccm

def mutual_information( list1, list2 ):
  contingency_matrix = clccm.contingency_matrix( list1, list2 )
  contingency = np.array( contingency_matrix, dtype = 'float' ) 
  contingency_sum = np.sum( contingency )
  pi = np.sum( contingency, axis = 1 )
  pj = np.sum( contingency, axis = 0 )
  outer = np.outer( pi, pj )
  nnz = contingency != 0.0
  contingency_nm = contingency[nnz]
  log_contingency_nm = np.log( contingency_nm )
  contingency_nm /= contingency_sum
  log_outer = -np.log( outer[nnz] ) + log( pi.sum() ) + log( pj.sum() )
  mi = ( contingency_nm * ( log_contingency_nm - log( contingency_sum ) )
        + contingency_nm * log_outer )
  mi = mi.sum()
  return mi
  
def randIndex( list1, list2 ):
  contingency_matrix = clccm.contingency_matrix( list1, list2 )
  Nsample = len(list1)
  sum_c1 = sum( comb( n_c, 2, exact = 1 ) for n_c in contingency_matrix.sum( axis = 1 ) )
  sum_c2 = sum( comb( n_k, 2, exact = 1 ) for n_k in contingency_matrix.sum( axis = 0 ) )

  sum_ = sum( comb( nij, 2, exact = 1 ) for nij in contingency_matrix.flatten() )
  prod_ = ( sum_c1 * sum_c2 ) / float( comb( Nsample, 2 ) )
  mean_ = ( sum_c2 + sum_c1 ) / 2.
  return ( ( sum_ - prod_ ) / ( mean_ - prod_ ) )
  
def jacard( list1, list2 ):
  pass