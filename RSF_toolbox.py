#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: 
"""

from toolbox import Truth_table_entry, bool_list_to_integer

"""*********************************************************************
*********************************Tool Box*******************************
*********************************************************************"""

def compute_representatives(locality):
    """
    Returns the representatives for a given locality.
    The output contains:\n
        -the representatives sorted by weight, then by lexicographic order.\n
        -a list of indexes to find the representative of a non-representative element.

    Parameters
    ----------
    locality : integer
        number of variables to consider.

    Returns
    -------
    (representatives, truth_table_index) : (array of arrays of Booleans, array of integer)
    The output contains:\n
        -the representatives sorted by weight, then by lexicographic order.\n
        -a list of indexes to find the representative of a non-representative element.
        For example, if truth_table_index[6] = 3, then the representative of the 8th element is the 3rd (starting at zero).

    """
    
    inputs = Truth_table_entry(locality)
    representatives = []
    truth_table_index = [-1]*(2**locality)
    
    for i in range(2**locality):   #for all entries
        new_input = [j for j in inputs.next_entry()] #local copy
        is_representative = True
        for rotate in range(locality-1): #all possible rotations
            new_input.insert(0, new_input.pop())    #rotate
            if bool_list_to_integer(new_input) < i: #it is not a representative
                is_representative = False
                truth_table_index[i] = bool_list_to_integer(new_input) #index to the representative
        if is_representative:
            new_input.insert(0, new_input.pop())    #rotate
            representatives.append(new_input)
    representatives.sort(key = sum) #sort by weight
    return (representatives, truth_table_index)

def representative_to_ANF(representative):
    """
    Compute the list of monomials corresponding to a representative.

    Parameters
    ----------
    representative : array of Booleans
        the representative to convert.

    Returns
    -------
    string : string of characters
        printable ANF.

    """
    l = len(representative)
    index_1 = []
    list_monomials = []
    
    
    #indexes of 1 in the representative
    for i in range(l):
        if representative[i] == 1:
            index_1.append(i)
            
    #building all monomials by rotation
    stop = False
    offset = 0        
    while not stop:
        new_monomial = [(index + offset)%l for index in index_1]
        new_monomial.sort()
        if new_monomial == index_1 and offset>0: #all rotations covered
            stop = True
        else:
            list_monomials.append(new_monomial)
        offset += 1
                
    list_monomials.sort()
    
    #convert into string
    string = ""
    for monomial in list_monomials:
        if string != "":
            string += " + "
        for variable in monomial:
            string += "X" + str(variable) 
            
    if string == "":
        string = "1"
    
    return string
    

def vector_include(v1, v2):
    """
    Checks whether v1 includes v2. For example:\n
    -vector_include([1,1,1,0], [1,0,1,0]) -> True
    -vector_include([1,1,0,0], [1,0,1,0]) -> False

    Parameters
    ----------
    v1 : array of Booleans
        first vector.
    v2 : array of Booleans
        second vector.

    Returns
    -------
    Boolean
        True if v1 includes v2. False otherwise.

    """
    
    
    for i in range(len(v1)):
        if v1[i]==0 and v2[i]==1:
            return False
    return True

def vector_orbit(v):
    """
    Compute the orbit of a vector.

    Parameters
    ----------
    v : array of Booleans
        vector to rotate.

    Returns
    -------
    orbit : array of Boolean vectors
        The orbit of the vector, obtained by rotation.

    """
    orbit = []
    new_vector = [i for i in v]
    vector_to_push = [i for i in v]
    orbit.append(vector_to_push)
    
    new_vector.insert(0, new_vector.pop())
    while new_vector != v:
        vector_to_push = [i for i in new_vector]
        orbit.append(vector_to_push)
        new_vector.insert(0, new_vector.pop())
    
    return orbit
    
def build_SANF_to_STT(representatives):
    """
    Conversion matrix from SANF to STT.
    """
    n = len(representatives)
    SANF_to_STT = [[0]*n for i in range(n)] #null matrix in F2^(n*n)
    
    for i in range(n): #for each row, representatives[i]
        orbit = vector_orbit(representatives[i])
        for j in range(n): #for each column, representatives[j]
            for vector in orbit: #for each vector in the orbit of representatives[i]
                if vector_include(representatives[j], vector):
                    SANF_to_STT[i][j] ^= 1
                    
    return SANF_to_STT

def build_STT_to_SWS(representatives):
    """
    Conversion matrix from STT to SWS.
    """
    n = len(representatives)
    STT_to_SWS = [[0]*n for i in range(n)] #null matrix in F2^(n*n)
    
    for i in range(n):  #for each row, representatives[i]
        orbit = vector_orbit(representatives[i])
        for j in range(n): #for each column, representatives[j]
            for vector in orbit:    #for each vector in the orbit of representatives[i]
                STT_to_SWS[i][j] += (-1)**(sum([representatives[j][k] * vector[k] for k in range(len(vector))])) #(-1) ** (scalar)
                
    return STT_to_SWS




coverage_locality = 0
representative_coverage = []

def representative_cover(locality, representative_index1, representative_index2):
    """
    Checks whether representative indexed by representative_index1 covers an element
    of the orbit of the representative indexed by representative_index2.\n
    The first call to this function creates a lookup table, and therefore can be slow.

    Parameters
    ----------
    locality : integer
        number of variables.
    representative_index1 : integer
        first representative index.
    representative_index2 : integer
        second representative index.

    Returns
    -------
    Boolean
        True if the first representative covers the second.

    """
    global coverage_locality
    global representative_coverage
    
    #initialize the representative_coverage matrice if necessary
    if coverage_locality != locality:
        (representatives, dummy) = compute_representatives(locality)
        nb_representatives = len(representatives)
        representative_coverage = [[False]*nb_representatives for i in range(nb_representatives)]
        for i in range(nb_representatives):
            for j in range(nb_representatives):
                if i == j:
                    representative_coverage[i][j] = True
                elif sum(representatives[i]) <= sum(representatives[j]): #if HW of j is greater than HW of i, i cannot cover j.
                    representative_coverage[i][j] = False
                else:
                    #default is False
                    for v in vector_orbit(representatives[j]): #for every vector in the orbit of j
                        if vector_include(representatives[i],v): #if i covers the vector
                            representative_coverage[i][j] = True 
        coverage_locality = locality
    #end initialize
    
    return representative_coverage[representative_index1][representative_index2]
    

def representative_distance(representative1, representative2):
    """
    Computes the minimal Hamming distance between two orbits.

    Parameters
    ----------
    representative1 : array of Booleans
        first orbit specified by a representative.
    representative2 : array of Booleans
        second orbit specified by a representative.

    Returns
    -------
    distance_min : integer
        minimal Hamming distance.

    """
    distance_min = len(representative1)
    for r in vector_orbit(representative1):
        tmp = [r[i]^representative2[i] for i in range(len(r))]
        distance_min = min(distance_min, sum(tmp))
    return distance_min












