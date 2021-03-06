#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 16:15:43 2019

@author: 
"""

"""*********************************************************************
*********************************Tool Box*******************************
*********************************************************************"""

from copy import deepcopy, copy

class Truth_table_entry:
    """
    Class Truth_table_entry.\n
    It is used as an iterator for a truth table, an ANF or a Walsh Spectrum.\n
    Consecutive calls to the next_entry() method with return: [0,0,0,0] -> [0,0,0,1] -> [0,0,1,0]...
    """
    
    current = []
    locality = 0
    def __init__(self, locality):
        """
        Constructor

        Parameters
        ----------
        locality : integer
            number of variables.

        Returns
        -------
        None.

        """
        self.current = [1]*locality
        self.locality = locality
    
    def next_entry(self):
        """
        Computes and returns the next entry.

        Returns
        -------
        array of Booleans
            the next entry of the truth table.

        """
        i = self.locality-1
        while self.current[i] == 1 and i >= 0:
            self.current[i] = 0
            i -= 1
        if i >= 0:
            self.current[i] = 1
        return self.current
    
    def previous_entry(self):
        """
        Computes and returns the previous entry.

        Returns
        -------
        array of Booleans
            the previous entry of the truth table.

        """
        i = self.locality-1
        while self.current[i] == 0 and i >= 0:
            self.current[i] = 1
            i -= 1
        if i >= 0:
            self.current[i] = 0
        return self.current

def bool_list_to_integer(bool_list):
    """
    Converts an array of Booleans into an integer.

    Parameters
    ----------
    bool_list : array of Booleans
        DESCRIPTION.

    Returns
    -------
    integer

    """
    ret = 0
    for i in range(len(bool_list)):
        ret = (ret<<1) + bool_list[i]
    return ret


def rank_increase(matrix):
    """
    Checks whether the last column of the given Boolean matrix increases the rank.\n 
    The matrix must be given as an array of columns.\n 
    The last column of the given matrix is modified by this function.\n 
    This function MUST be called for each column added to the matrix, or the result might be incorrect.

    Parameters
    ----------
    matrix : array of arrays of Booleans
        the Boolean matrix to check, seen as an array of columns.

    Returns
    -------
    bool
        True if the last columns increases the rank by one, False otherwise.

    """
    nb_columns = len(matrix)
    nb_rows = len(matrix[0])
    last_column = nb_columns - 1
    
    for i in range(nb_columns-1): #for each non-last column
        position_premier_1 = 0
        while position_premier_1 < nb_rows and matrix[i][position_premier_1] == 0: #research the first 1 of the column
            position_premier_1 += 1

        if position_premier_1 < nb_rows and matrix[last_column][position_premier_1] == 1:   #if the last column has the same 1, cancel it.
            for j in range(position_premier_1, nb_rows):
                matrix[last_column][j] = matrix[last_column][j] ^ matrix[i][j]
    
    #if the last column is null, the rank is not increased
    if matrix[last_column] == [0 for i in range(nb_rows)]:
        return False 
    else:
        return True


def moebius_poly(moebius_table,l):
    """
    Convert the ANF as a printable string.
    """
    inputs = Truth_table_entry()
    out = 0
    for i in range(len(moebius_table)):
        monome = ""
        inp = inputs.next_entry()
        if moebius_table[i]==1:
            for k in range(0, l):
                if inp[k] == 1 :
                    monome = monome+"X" + str(l-1-k)
            if out == 0:
                if monome == "":
                    out = "1"
                else:
                    out = monome
            else:
                out = monome + " + " + out
    out = "P = " + out
    return out

def sign(f):
    """
    Function sign. 0 is mapped to -1, and 1 is mapped to 1.

    Parameters
    ----------
    f : array of Booleans
        truth table of a functions.

    Returns
    -------
    array of {-1;1}
        sign applied to each element of the truth table.

    """
    res = []
    for e in f:
        res.append(2*e-1)
    return res

def Walsh_transform(f,locality):
    """
    Walsh transform of a function.

    Parameters
    ----------
    f : TYPE
        DESCRIPTION.
    locality : TYPE
        DESCRIPTION.

    Returns
    -------
    F : TYPE
        DESCRIPTION.

    """
    
    F = copy(f)
    divide = (2**locality)>>1
    while divide != 0:
        offset = divide
        for conquer in range(0,(2**locality),divide<<1):
            for entry in range(conquer, conquer+offset):
                tmp = F[entry]
                F[entry] = F[entry] + F[entry + offset]
                F[entry + offset] = tmp - F[entry + offset]
        divide = divide >>1
    return F

def Moebius_transform(f,locality):
    """
    Moebius transform to convert ANF into truth table and vice-versa.

    Parameters
    ----------
    f : array of Booleans
        truth table (or equivalently ANF).
    locality : integer
        locality

    Returns
    -------
    F : array of Booleans
        ANF (or equivalently truth table).

    """
    F = copy(f)
    divide = (2**locality)>>1
    while divide != 0:
        offset = divide
        for conquer in range(0,(2**locality),divide<<1):
            for entry in range(conquer, conquer+offset):
                F[entry + offset] = F[entry] ^ F[entry + offset]
        divide = divide >>1
    return F

def truncate_degree(f, deg_min, deg_max):
    """
    Truncates the specified degree (resp. weight) of a specified ANF (resp. Walsh spectrum).

    Parameters
    ----------
    f : array of Booleans or integers
        ANF or Walsh spectrum.
    deg_min : integer
        minimal degree or weight.
    deg_max : integer
        maximal degree or weight.

    Returns
    -------
    res : array of Booleans or integers
        truncated result.

    """
    res = []
    for i in range(len(f)):
        deg = 0
        j=0
        while (i>>j) != 0:
            if (i>>j)%2 == 1:
                deg += 1
            j += 1
        if deg >= deg_min and deg <= deg_max:
            res.append(f[i])
    return res

def vector_complement(v):
    vp = [0]*len(v)
    for i in range(len(v)):
        if v[i] == 0:
            vp[i] = 1
    return vp

def vector_opposite(v):
    vp = [0]*len(v)
    for i in range(len(v)):
        vp[i] = -v[i]
    return vp

def xor_matrices(matrix1, matrix2):
    nb_rows = len(matrix1)
    nb_columns = len(matrix1[0])
    sum_matrix = deepcopy(matrix1)
    
    for i in range(nb_rows):
        for j in range(nb_columns):
            sum_matrix[i][j] ^= matrix2[i][j]
    return sum_matrix

def add_vectors(v1, v2):
    sum_v = copy(v1)
    
    for i in range(len(v1)):
        sum_v[i] += v2[i]
    return sum_v

def xor_vectors(v1, v2):
    sum_v = copy(v1)
    
    for i in range(len(v1)):
        sum_v[i] ^= v2[i]
    return sum_v

def permute(x, permutation):
    """
    Permutes the variables of x following the given permutation.
    The permutation is an array of indexes of length the number of variables of x.\n
    Example:\n
    permutation(9,[0,3,1,2]) outputs 12 after applying the following permutation:\n
    Y0 <- X0\n
    Y1 <- X3\n
    Y2 <- X1\n
    Y3 <- X2\n

    Parameters
    ----------
    x : integer
    permutation : array of integers

    Returns
    -------
    y : integer

    """
    nb_variables = len(permutation)
    y = 0
    for i in range(len(permutation)):
        y = y << 1
        y |= (x & (1 << (nb_variables-1-permutation[i]))) and 1
    return y

        