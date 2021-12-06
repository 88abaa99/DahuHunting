#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: 
    
This module defines a BF (Boolean Function) class.
An instantiation represents a Boolean function.

This module allows to work alternatively on its truth table, Walsh spectrum or ANF,
while keeping the two other representations updated.

This module also allows to check its resiliency and algebraic immunity,
and to compute its annihilators.
"""

from copy import deepcopy as copy
from toolbox import sign, Walsh_transform, Moebius_transform, vector_complement, permute
from annihilator import get_annihilators

"""*********************************************************************
    Global variables
*********************************************************************"""



"""*********************************************************************
    Class BF
*********************************************************************"""
class BF:
    l = 0   #locality
    
    WS = []
    is_WS_uptodate = False
    
    TT = []
    is_TT_uptodate = False
    
    ANF = ""
    is_ANF_uptodate = False
    
    monomials_degree = []   #degree/weight of each monomial/element for the ANF and WS
    
    annihilators_basis_f = []   #annihilators of f
    annihilators_basis_fp1 = []   #annihilators of f+1

    
    
    
    """*********************************************************************
    Constructor
    *********************************************************************"""
    def __init__(self, locality):
        """
        Constructor of the class BF (Boolean Function).

        Parameters
        ----------
        locality : integer
            Locality of the Boolean function

        Returns
        -------
        None.

        """
        
        self.l = locality
        
        self.WS = [0]*(2**self.l)
        self.TT = [0]*(2**self.l)
        self.ANF= [0]*(2**self.l)
        
        self.is_WS_uptodate = False
        self.is_ANF_uptodate = True
        self.is_TT_uptodate = True
        
        #degree of the monomials
        self.monomials_degree = []
        for x in range(2**self.l):
            deg = 0
            while x != 0:
                if x%2 :
                    deg += 1
                x >>= 1
            self.monomials_degree.append(deg)
            
        #basis of annihilators
        self.annihilators_basis_f = []
        self.annihilators_basis_fp1 = []
        
        return
        
    def set_ANF(self, new_SANF):
        """
        Set a new ANF.

        Parameters
        ----------
        new_TT : Array of Booleans
            The new ANF of the function.

        Returns
        -------
        None.

        """
        self.ANF = copy(new_SANF)
        self.is_WS_uptodate = False
        self.is_ANF_uptodate = True
        self.is_TT_uptodate = False
        return
    
    def set_TT(self, new_TT):
        """
        Set a new truth table.

        Parameters
        ----------
        new_TT : Array of Booleans
            The new truth table of the function.

        Returns
        -------
        None.

        """
        self.TT = copy(new_TT)
        self.is_WS_uptodate = False
        self.is_ANF_uptodate = False
        self.is_TT_uptodate = True
        return

    def update_TT_from_ANF(self):
        """
        Update the truth table with the ANF.
        
        Requirements
        ------------
        The ANF must be up to date.

        Returns
        -------
        None.

        """
        self.TT = Moebius_transform(self.ANF, self.l)
        return
    
    def update_WS_from_TT(self):
        """
        Update the Walsh spectrum with the truth table.
        
        Requirements
        ------------
        The truth table must be up to date.

        Returns
        -------
        None.

        """
        self.WS = Walsh_transform(sign(self.TT), self.l)
        return
    
    def update_ANF_from_TT(self):
        """
        Update the ANF with the truth table.
        
        Requirements
        ------------
        The truth table must be up to date.

        Returns
        -------
        None.

        """
        self.ANF = Moebius_transform(self.TT, self.l)
        return
    

    def update_TT(self):
        """
        Update the truth table.
        
        Future developpements
        ---------------------
        So far, it can only be updated from the ANF. Update from WS will be added soon.

        Returns
        -------
        None.

        """
        if self.is_TT_uptodate: #if TT is up to date, nothing to do
            return
        if self.is_ANF_uptodate:  #update with the ANF
            self.update_TT_from_ANF()
            self.is_TT_uptodate = True
            return
        return
    

    def update_WS(self):
        """
        Update the Walsh Spectrum.
        
        Future developpements
        ---------------------
        So far, it can only be updated from the truth table. Update from ANF will be added soon.

        Returns
        -------
        None.

        """
        if self.is_WS_uptodate: #if WS is up to date, nothing to do
            return

        self.update_TT() #first, update TT
        
        self.update_WS_from_TT() #then update with TT
        self.is_WS_uptodate = True
        return
    
    def update_ANF(self):
        """
        Update the ANF.
        
        Future developpements
        ---------------------
        So far, it can only be updated from the truth table. Update from WS will be added soon.

        Returns
        -------
        None.

        """
        if self.is_ANF_uptodate: #if ANF is up to date, nothing to do
            return
        if self.is_TT_uptodate:  #update with the TT
            self.update_ANF_from_TT()
            self.is_ANF_uptodate = True
            return
        return

    def string_ANF(self):
        """
        Returns the ANF as a printable string.
        
        Requirements
        ------------
        The ANF must be up to date.

        Returns
        -------
        string : string of characters
            ANF as a printable string of characters.

        """
        string = ""
        if self.ANF[0] == 1:
            string = "1"
            
        for monome in range(1,2**self.l):
            if self.ANF[monome] == 1:
                if string != "":
                    string += " + "
                index = 0
                while monome >> index != 0:
                    if (monome>>index)%2 == 1:
                        string += "X"+str(index)
                    index += 1
        return string
    
    def update_annihilators(self, max_degree):
        """
        Update the basis of annihilators of f and f+1 for a given maximal degree.
        
        Requirements
        ------------
        The Walsh spectrum must be up to date.

        Parameters
        ----------
        max_degree : integer
            Maximal degree of the annihilators

        Returns
        -------
        None.

        """
        self.annihilators_basis_f = [] #delete the old basis of annihilators of f
        basis = get_annihilators(self.TT, self.l, max_degree) #compute the new basis in the form [int,int,...]
        
        for annihilator in basis: #for each annihilator of this basis, in the form [int, int, ...]
            annihilator_ANF = [0]*(2**self.l)
            for index in annihilator: #convert it in ANF [0,0,1,...]
                annihilator_ANF[index] = 1
            self.annihilators_basis_f.append(annihilator_ANF)
            
        self.annihilators_basis_fp1 = [] #delete the old basis of annihilators of f1
        basis = get_annihilators(vector_complement(self.TT), self.l, max_degree) #compute the new basis in the form [int,int,...]
        
        for annihilator in basis: #for each annihilator of this basis, in the form [int, int, ...]
            annihilator_ANF = [0]*(2**self.l)
            for index in annihilator: #convert it in ANF [0,0,1,...]
                annihilator_ANF[index] = 1
            self.annihilators_basis_fp1.append(annihilator_ANF)
        return
    
    def get_annihilators(self, min_degree=0, max_degree=-1):
        """
            Returns the basis of annihilators of f and f+1 in their ANF form.\n
            If min_degree is specified, only monomials of degree greater or equal to min_degree are returned.\n
            If max_degree is specified, only monomials of degree smaller or equal to max_degree are returned.
            
        Parameters
        ----------
        min_degree : integer, optional
            Minimal degree to truncate annihilators.
        max_degree : integer, optional
            Maximal degree to truncate annihilators.

        Returns
        -------
        (Array, Array)
            Two arrays of annihilators in their ANF representation.

        """
        if min_degree==0 and max_degree==-1: #no degree truncature
            return copy(self.annihilators_basis)
        
        #annihilators of f
        annihilator_list_f = []
        for annihilator in self.annihilators_basis_f: #for each annihilator
            new_annihilator = []
            for index in range(2**self.l): #keep only the correct degrees
                if self.monomials_degree[index] >= min_degree and self.monomials_degree[index] <= max_degree:
                    new_annihilator.append(annihilator[index])
            annihilator_list_f.append(new_annihilator)
            
        #annihilators of f+1
        annihilator_list_fp1 = []
        for annihilator in self.annihilators_basis_fp1: #for each annihilator
            new_annihilator = []
            for index in range(2**self.l): #keep only the correct degrees
                if self.monomials_degree[index] >= min_degree and self.monomials_degree[index] <= max_degree:
                    new_annihilator.append(annihilator[index])
            annihilator_list_fp1.append(new_annihilator)
            
        return (annihilator_list_f, annihilator_list_fp1)
    
    def get_WalshSpectrum(self, min_weight=0, max_weight=-1):
        """
            Returns the Walsh spectrum of the function.\n
            If min_weight is specified, only elements of weight greater or equal to min_weight are returned.\n
            If max_weight is specified, only elements of weight smaller or equal to max_weight are returned.
        
        Parameters
        ----------
        min_weight : integer, optional
            Minimal weight of the returned elements.
        max_weight : integer, optional
            Maximal weight of the returned elements.

        Returns
        -------
        Array
            Walsh spectrum as an array of integer.
            
        """
        if min_weight==0 and max_weight==-1: #no weight truncature
            return copy(self.WS)
        
        new_WS = []
        for index in range(2**self.l): #keep only the asked weights
            if self.monomials_degree[index] >= min_weight and self.monomials_degree[index] <= max_weight:
                new_WS.append(self.WS[index])
        return new_WS
    
    def get_ANF(self, min_degree=0, max_degree=-1):
        """
            Returns the ANF of the function.\n
            If min_degree is specified, only monomials of degree greater or equal to min_degree are returned.\n
            If max_degree is specified, only monomials of degree smaller or equal to max_degree are returned.
        
        Parameters
        ----------
        min_degree : integer, optional
            Minimal degree to truncate the ANF.
        max_degree : integer, optional
            Maximal degree to truncate the ANF.

        Returns
        -------
        Array
            ANF as an array of Booleans.
            
        """
        if min_degree==0 and max_degree==-1: #ne degree truncature
            return copy(self.ANF)
        
        new_ANF = []
        for index in range(2**self.l): #keep only the correct degrees
            if self.monomials_degree[index] >= min_degree and self.monomials_degree[index] <= max_degree:
                new_ANF.append(self.ANF[index])
        return new_ANF
    
    def is_resilient(self, r):
        """
        Returns True if the Boolean function is r-resilient, False otherwise.\n
        Warning: the Walsh spectrum must be up-to-date.

        Parameters
        ----------
        r : integer
            resilience to verify.

        Returns
        -------
        bool
            True if the Boolean function is r-resilient.

        """
        for i in range(len(self.WS)):
            if (self.monomials_degree[i] <= r) and (self.WS[i] != 0):
                return False
        return True
    
    def is_algebraic_immune(self, ai):
        """
        Returns True if the Boolean function is ai-algebraic-immune, False otherwise.\n
        Warning: the truth table spectrum must be up-to-date.

        Parameters
        ----------
        ai : integer
            algebraic immunity to verify.

        Returns
        -------
        bool
            True if the Boolean function is ai-algebraic-immune.

        """
        #annulateurs de f+1
        if get_annihilators(self.TT, self.l, ai-1) != []:
            return False
        
        #annulateurs de f+1
        if get_annihilators(vector_complement(self.TT), self.l, ai-1) != []:
            return False
        
        return True
    
    def permute(self, permutation, new_object = False):
        """
        Permutes the variables of the function following the given permutation.
        The permutation is an array of indexes of length the number of variables of x.\n
        Example:\n
        for a function f of locality 4, f.permute([0,3,1,2]) gives:\n
        f'(X0,X1,X2,X3) = f(X0,X3,X1,X2)
        
        If new_object is set to False (default), then the calling object is directly modified.\n
        If new_object is set to True, the calling object is not modified and a new object is created and returned.

        Parameters
        ----------
        permutation : array of integers

        Returns
        -------
        None if new_object is set to False (default),\n
        a new object with the permuted function otherwise. 

        """
        new_TT = [0]*(2**self.l)
        for i in range(2**self.l):
            new_TT[i] = self.TT[permute(i, permutation)]
            
        if new_object == False:
            self.TT = new_TT
            self.is_TT_uptodate = True
            self.is_ANF_uptodate = False
            self.is_WS_uptodate = False
            return
        else:
            new_bf = BF(self.l)
            new_bf.set_TT(new_TT)
            return new_bf
        
    def translate(self, translation, new_object = False):
        """
        Translate the Boolean function according to the given translation.\n
        A translation vector t maps a function f(x) to f(x^t).\n
        If new_object is set to False (default), then the calling object is directly modified.\n
        If new_object is set to True, the calling object is not modified and a new object is created and returned.

        Parameters
        ----------
        translation : Boolean vector of l elements.

        Returns
        -------
        None if new_object is False (default).\n
        Otherwise, a new instance of BF is returned.

        """
        new_TT = [0]*(2**self.l)
        mask = 0
        for t in translation:
            mask = (mask << 1) ^ t
        for i in range(2**self.l):
            new_TT[i] = self.TT[mask ^ i]
        if new_object == False:
            self.TT = new_TT
            self.is_TT_uptodate = True
            self.is_ANF_uptodate = False
            self.is_WS_uptodate = False
            return
        else:
            new_bf = BF(self.l)
            new_bf.set_TT(new_TT)
            return new_bf