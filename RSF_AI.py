#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 2020

This module contains the class Verification_AI,
that allows to check the algebraic immunity of function defined as a truth table.

@author:
"""

from reedmuller import ReedMuller
from toolbox import rank_increase, bool_list_to_integer

reedmuller_param = (0,0)
RM = []

def get_reedmuller(r,m):
    """
    Returns the Reed-Muller matrix RM(r,m)

    Parameters
    ----------
    r : integer
        order.
    m : integer
        number of variables.

    Returns
    -------
    RM : matrix of Booleans
        RM(r,m).

    """
    global reedmuller_param
    global RM
    
    #if not the same parameters as last call
    if reedmuller_param != (r,m):
        RM = ReedMuller(r,m)
        reedmuller_param = (r,m)
    return RM

class Verification_AI:
    """
    The Verification_AI class allows to check the algebraic immunity of a function.\n
    The function must be defined as a truth table and the check_and_add method must be called for each element of the truth table.\n
    
    """
    
    l = 0
    ai = 0
    M_AI = []
    Mat = [[],[]]
    nb_remaining_elements = 0 #number of elements to add to the truth table
    rank = [0,0]
    rank_max = 0
    nb_monomes_AI = 0

    def __init__(self, locality, algebraic_immunity):
        """
        Constructor
        
        Parameters
        ----------
        locality : integer
            number of variables to consider.
        algebraic_immunity : integer
            algebraic immunity to check.
        
        Returns
        -------
        None.
        """
        
        self.l = locality
        self.ai = algebraic_immunity
        self.Mat = [[],[]]
        self.rank = [0,0]
        self.nb_remaining_elements = 2**(self.l)
        
        #Initialisation of RM(r,m)
        M_AI = get_reedmuller(algebraic_immunity-1,locality)
        M_AI_tmp = [[M_AI[i][j] for j in range(len(M_AI[i]))] for i in range(len(M_AI))]    #conversion
        self.M_AI = M_AI_tmp
        
        self.nb_monomes_AI = len(M_AI[0])
        self.rank_max = self.nb_monomes_AI
        
    def corresponding_line_AI(self, l_uple):
        i = 2**self.l - (1 + bool_list_to_integer(l_uple))
        return self.M_AI[i][0:self.nb_monomes_AI]

    def check_and_add(self, X, y):   
        """
        Checks whether adding f(X) = y to the current truth table violates the specified algebraic immunity.\n
        If not, the method returns True, otherwise it returns False.\n
        To verify an entire an entire truth table, this method must be called for each element (in any order).

        Parameters
        ----------
        X : array of Booleans
            input X of the function.
        y : Boolean
            evaluation of f in X.

        Returns
        -------
        Boolean
            True if the function can still be reach the specified algebraic immunity.\n
            False otherwise.

        """
        self.nb_remaining_elements -= 1
        
        #recover the corresponding line in RM(r,m)
        line = self.corresponding_line_AI(X)
        
        #add it to the matrix
        self.Mat[y].append(line)
        
        #check whether it increases the rank
        res = rank_increase(self.Mat[y])
        if res == True: #rank default
            self.rank[y] += 1
        
        if sum(self.rank) + self.nb_remaining_elements < 2*self.rank_max: #there exists an annihilator of lower degree
            return False
        else:
            return True

    def reset(self):
        """
        Reset the current truth table to null.\n
        This method should be preferred to re-instantiating the class if the algebraic immunity target is the same.

        Returns
        -------
        None.

        """
        self.Mat = [[],[]]
        self.rank = [0,0]
        self.nb_remaining_elements = 2**(self.l)
        return
