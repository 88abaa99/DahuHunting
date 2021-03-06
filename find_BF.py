#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: 
"""

from BF import BF
from toolbox import Truth_table_entry


def find_BF_naive(locality, resiliency, algebraic_immunity):
    """
    Naive approach to find Boolean functions of a specified locality, resiliency and algebraic immunity.
    The given resiliency and algebraic immunity are not strict but lower bounds.

    Parameters
    ----------
    locality : integer
        number of variables to consider.
    resilience : integer
        minimal resiliency.
    algebraic_immunity : integer
        minimal algebraic immunity.

    Returns
    -------
    found : integer
        number of functions satisfying the criteria.

    """
    bf = BF(locality)
    tt = Truth_table_entry(2**locality)
    found = 0
    for i in range(2**(2**locality)):
        bf.set_TT(tt.next_entry())
        bf.update_WS()
        if (resiliency==-1 or bf.is_resilient(resiliency)) and bf.is_algebraic_immune(algebraic_immunity):
            bf.update_TT()
            found += 1
            print(bf.TT)
       
    return found