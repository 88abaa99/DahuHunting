

"""
Created on Thu Mar 28 15:00:00 2019

This module allows to build a Reed-Muller matrix.

@author:
"""

import itertools
from toolbox import Truth_table_entry

def ReedMuller(r,m):
    n = 2**m
    monomials = []
    RM = []
    
    #Initialising monomials
    for i in range(r+1):
        for monomial in itertools.combinations(range(m),i):
            monomials.append(monomial)
                
    #Row computation
    entry = Truth_table_entry(m)
    for i in range(n):
        row = []
        for monomial in monomials:
            check = 1
            for variable in monomial:
                if entry.current[variable] == 0:
                    check = 0
            row.append(check)
        entry.previous_entry()
        RM.append(row)
    
    return RM
