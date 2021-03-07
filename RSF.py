#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author:

This module defines an RSF (Rotational Symmetric Function) class.
An instantiation represents a rotationnal symmetric function.

This module allows to work alternatively on its STT, SWS or SANF.
while keeping the two other representations updated.

Limitation : so far, we can only work on the SANF, the two others are updated, but we can not modify them.

This module also allows to check its resiliency and algebraic immunity.
"""

from toolbox import Truth_table_entry, bool_list_to_integer
from RSF_toolbox import compute_representatives, representative_to_ANF, build_SANF_to_STT, build_STT_to_SWS
from RSF_AI import Verification_AI


class RSF:
    l = 0   #locality
    
    SANF = []
    is_SANF_uptodate = False
    
    STT = []
    is_STT_uptodate = False
    SANF_to_STT = []
    
    SWS = []
    is_SWS_uptodate = False
    build_STT_to_SWS = []
    
    TT = []
    is_TT_uptodate = False
    STT_to_TT = []
    
    ANF = ""
    is_ANF_uptodate = False
    SANF_to_ANF = []
    
    
    representatives = []    #representatives sorted by weight, then by increasing order
    nb_representatives = 0  #number of representatives
    nb_representatives_by_weight = []   #number of representatives by weight
    
    verification_AI = 0
    
    

    def __init__(self, locality):
        """
        Constructor of the class RSF (Rotational Symmetric Function).

        Parameters
        ----------
        locality : integer
            Locality of the rotational symmetric function

        Returns
        -------
        None.

        """
        self.l = locality
        (self.representatives, self.STT_to_TT) = compute_representatives(locality)
        self.nb_representatives = len(self.representatives)
        
        
        self.SANF = [0]*self.nb_representatives
        self.STT = [0]*self.nb_representatives
        self.SWS = [0]*self.nb_representatives
        self.TT = [0]*(2**self.l)
        self.ANF = "0"
        
        self.is_SANF_uptodate = True
        self.is_STT_uptodate = True
        self.is_SWS_uptodate = False
        self.is_ANF_uptodate = True
        self.is_TT_uptodate = True
        
        self.SANF_to_ANF = ['('+representative_to_ANF(r)+')' for r in self.representatives]
        self.SANF_to_STT = build_SANF_to_STT(self.representatives)
        self.STT_to_SWS = build_STT_to_SWS(self.representatives)
        
        self.verification_AI = Verification_AI(locality, int((locality+1)/2)) #optimal AI by default
        
        self.nb_representatives_by_weight = [0]*(locality+1)
        for r in self.representatives:
            self.nb_representatives_by_weight[sum(r)] += 1
        
        return
        
        
    def set_SANF(self, new_SANF):
        """
        Set a new SANF (Simplified Algebraic Normal Form).

        Parameters
        ----------
        new_SANF : Array of Booleans
            The new SANF of the rotational symmetric function.

        Returns
        -------
        None.

        """
        self.SANF = new_SANF
        self.is_SANF_uptodate = True
        self.is_STT_uptodate = False
        self.is_SWS_uptodate = False
        self.is_ANF_uptodate = False
        self.is_TT_uptodate = False
        return

            
    def update_STT_from_SANF(self):
        """
        Update the STT with the SANF.

        Returns
        -------
        None.

        """
        
        #update STT
        self.STT = [0]*self.nb_representatives
        
        #multiply the vector SANF with the conversion matrix
        for i in range(self.nb_representatives):
            if self.SANF[i] == 1:
                for j in range(self.nb_representatives):
                    self.STT[j] ^= self.SANF_to_STT[i][j]
        
        self.is_STT_uptodate = self.is_SANF_uptodate
        return
    
    def update_SWS_from_STT(self):
        """
        Update the SWS with the STT.

        Returns
        -------
        None.

        """
        
        #initialisation
        self.SWS = [0]*self.nb_representatives
        
        #calcul de 1-2STT
        tmp = [1 - 2*self.STT[i] for i in range(self.nb_representatives)]
        
        #multiply the vector tmp with the conversion matrix
        for i in range(self.nb_representatives):
            for j in range(self.nb_representatives):
                self.SWS[j] += tmp[i] * self.STT_to_SWS[i][j]
        
        self.is_SWS_uptodate = self.is_STT_uptodate
        return
    

    def update_ANF_from_SANF(self):
        """
        Update the ANF with the SANF.

        Returns
        -------
        None.

        """
        #initialisation
        self.ANF = ""
        
        for i in range(self.nb_representatives):
            if self.SANF[i] == 1:
                if self.ANF != "":
                    self.ANF += " + "
                self.ANF += self.SANF_to_ANF[i]
                
        self.is_ANF_uptodate = self.is_SANF_uptodate
        return
    
    def update_TT_from_STT(self):
        """
        Update the TT with the STT.

        Returns
        -------
        None.

        """
        #initialisation
        self.TT = [0]*(2**self.l)
        
        #copy the representatives
        for i in range(self.nb_representatives):
            self.TT[bool_list_to_integer(self.representatives[i])] = self.STT[i]
            
            
        #copy the other elements
        for i in range(2**self.l):
            if self.STT_to_TT[i] != -1: #if it is not a representative
                self.TT[i] = self.TT[self.STT_to_TT[i]] #same value as its representative
                
        self.is_TT_uptodate = self.is_STT_uptodate
        return


    def update_STT(self):
        """
        Update the STT (Simplified Truth Table).
        
        Future developpements
        ---------------------
        So far, it can only be updated from the SANF. Update from SWS will be added soon.

        Returns
        -------
        None.

        """
        if self.is_STT_uptodate: #STT it already updated
            return
        if self.is_SANF_uptodate:  #update with SANF
            self.update_STT_from_SANF()
            return
        return
    
    def update_SWS(self):
        """
        Update the SWS (Simplified Walsh Spectrum).

        Returns
        -------
        None.

        """
        if self.is_SWS_uptodate: #SWS it already updated
            return
        if not self.is_STT_uptodate:  #update with STT
            self.update_STT()
        self.update_SWS_from_STT()
        return
    
    def update_TT(self):
        """
        Update the truth table.

        Returns
        -------
        None.

        """
        if self.is_TT_uptodate: #TT it already updated
            return
        if not self.is_STT_uptodate:  #update with STT
            self.update_STT()
        self.update_TT_from_STT()
        return


    def is_resilient_optimised(self, resilience):
        """
        
        Returns True if the Boolean function is r-resilient, False otherwise.\n
        Otimised version of the method is_resilient. However, is_resilient should be preferred if the Walsh spectrum is already updated.

        Parameters
        ----------
        r : integer
            resilience to verify.

        Returns
        -------
        bool
            True if the Boolean function is r-resilient.

        """
        #update STT
        self.update_STT()
        
        #initialisation
        partial_SWS = [0]*self.nb_representatives
        
        #compute 1-2STT
        tmp = [1 - 2*self.STT[i] for i in range(self.nb_representatives)]
        
        #number of elements of SWS to compute
        stop = sum(self.nb_representatives_by_weight[0:(resilience+1)])
        
        
        for i in range(stop):
            for j in range(self.nb_representatives):
                partial_SWS[i] += tmp[j] * self.STT_to_SWS[j][i] #multiply vector tmp by the conversion matrix
            
            if partial_SWS[i] != 0:  #its Walsh element must be null
                return False
        
        return True
    

    def is_resilient(self, resiliency):
        """
        Returns True if the Boolean function is r-resilient, False otherwise.\n

        Parameters
        ----------
        r : integer
            resilience to verify.

        Returns
        -------
        bool
            True if the Boolean function is r-resilient.

        """
        
        #update Walsh spectrum
        self.update_SWS()
        
        i = 0
        while sum(self.representatives[i]) <= resiliency: #if the representative has a weight smaller or equal to the resiliency
            if self.SWS[i] != 0:    #its Walsh element must be null
                return False
            i += 1
        return True
    
    def is_algebraic_immune(self, algebraic_immunity):
        """
        Returns True if the Boolean function is ai-algebraic-immune, False otherwise.\n

        Parameters
        ----------
        ai : integer
            algebraic immunity to verify.

        Returns
        -------
        bool
            True if the Boolean function is ai-algebraic-immune.

        """
        
        if self.verification_AI.ai != algebraic_immunity: #if parameters are different from previous call
            self.verification_AI = Verification_AI(self.l, algebraic_immunity) #full initialisation
        else:
            self.verification_AI.reset()    #partial reinitialisation
        
        #update truth table
        self.update_TT()
        
        #add every entry one by one
        inputs = Truth_table_entry(self.l)
        for i in range(2**self.l):
            new_input = inputs.next_entry()
            if not self.verification_AI.check_and_add(new_input, self.TT[i]):    #if an entry contradicts the algebraic immunity
                return False
        return True