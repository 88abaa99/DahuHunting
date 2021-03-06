# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 10:29:57 2020

@author: 
    
    
This module allows to compute a basis of annihilators of a functions.

Unlike other modules of this projet, monomials are represented by integers instead of array or Booleans.
However, a simple binary decomposition makes the two representations equivalent.
"""

from RSF_toolbox import Truth_table_entry, bool_list_to_integer


def evaluate_function(f,x):
    """
    Evaluates in x a function described as a list of monomials, each monomial represented by an integer.

    Parameters
    ----------
    f : array of integers
        the function to evaluate.
    x : integer

    Returns
    -------
    res : Boolean
        f(x).

    """
    
    res = 0
    for monome in f:
        if (x & monome) == monome: #if the monomial evaulated in x gives 1
            res = res ^ 1
    return res

def clean_monomials(annihilator):
    """
    Searches and removes multiply defined monomials in an annihilator.

    Parameters
    ----------
    annihilator : array of integers
        annihilator, each integer represents a monomial.

    Returns
    -------
    The input is directly modified.

    """
    annihilator.sort()
    i=0
    while i+1 < len(annihilator):
        if annihilator[i] == annihilator[i+1]:
            annihilator.pop(i+1)
            annihilator.pop(i)
        else:
            i += 1


def get_annihilators(f, locality, deg_annihilator):
    """
    Computes a basis of annihilators of the function f if it exists.\n
    The degree of the returned annihilators are smaller or equal to deg_annihilator.\n
    Each annihilator is represented as an array of integer, each integer representing a monomial.

    Parameters
    ----------
    f : array of Booleans
        function to annihilate.
    locality : integer
        locality of f.
    deg_annihilator : integer
        maximal degree of the annihilators.

    Returns
    -------
    array of arrays of integers [[int,int,int],[int,int],...].
        basis of annihilators.

    """

    #list of monomials of weight less or equal to deg_annihilator
    tt = Truth_table_entry(locality)
    monomials = []
    for i in range(2**locality):
        if(sum(tt.next_entry()) <= deg_annihilator):
            monomials.append(bool_list_to_integer(tt.current)) 
    
    #initialisation
    S=[[0]]
    i=0
    if 1 in f:
        x=f.index(1)
    else:
        return [[monomials[j]] for j in range (len(monomials))]
    
    while 1:
        #empile?
        while (i+1) <len(monomials) and monomials[i+1] <= x:
            S.append([monomials[i+1]])
            i += 1
            
        #depile?
        S_index= len(S) - 1
        while S_index>=0  and evaluate_function(S[S_index],x) != 1:
            S_index-=1
        if S_index>=0:
            for j in range(S_index):
                if evaluate_function(S[j],x)==1:
                    S[j] += S[S_index]
                    clean_monomials(S[j])
            S.pop(S_index) 
        
        #jump?
        """if S == [] or S == [[0]]:
            i += 1
            if i >= len(monomials):
                return []
            S.append([monomials[i]])
            x = monomials[i]"""
            
        #loop?
        xp = x+1
        while xp < len(f) and f[xp] != 1:
            xp+=1
        if xp == len(f):
            x = xp
            break
        else:
            x = xp
    
    #empile?
    while (i+1) <len(monomials) and monomials[i+1] <= x:
        S.append([monomials[i+1]])
        i += 1
    
    return S

def annihilator_to_ANF(annihilator):
    """
    Converts an annihilators in array of integers format into a printable ANF.

    Parameters
    ----------
    annihilator : array of integers
        annihilator to convert.

    Returns
    -------
    ANF : string

    """
    ANF = ""
    for monomial in annihilator:
        if ANF != "":
            ANF += " + "
        i=0
        if monomial==0:
            ANF += "1"
        while monomial > 0:
            if monomial%2 == 1:
                ANF += "X"+str(i)
            monomial = monomial >> 1
            i += 1
    return ANF
    

def truncate_annihilator(annihilator, degree):
    """
    Truncates an annihilator to keep only a specified degree part.

    Parameters
    ----------
    annihilator : array of integers
        annihilator to truncate.
    degree : integer
        degree to keep.

    Returns
    -------
    truncated_annihilator : array of integers
        truncated annihilator.

    """
    
    truncated_annihilator = []
    for monome in annihilator:
        hw = 0
        i = 0
        while (monome >> i) != 0: #calcul du degré du monôme
            if (monome >> i) & 1 == 1:
                hw += 1
            i += 1
        if hw >= degree:
            truncated_annihilator.append(monome)
    truncated_annihilator.sort()
    return truncated_annihilator
        
        
        
        
        