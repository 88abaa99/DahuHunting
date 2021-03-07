#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author:
"""

from RSF import RSF
from RSF_toolbox import Truth_table_entry, bool_list_to_integer, representative_cover
import time, os

def find_RSF_from_SANF_naive(locality, resiliency, algebraic_immunity):
    """
    Naive exhaustive approach to find RSF with specified locality, resiliency and algebraic immunity.

    Parameters
    ----------
    locality : integer
        number of variables to consider.
    resiliency : integer
        minimal resiliency to verify.
    algebraic_immunity : integer
        minimal algebraic immunity to verify.

    Returns
    -------
    found : integer
        number of results.

    """
    rsf = RSF(locality)
    tt = Truth_table_entry(rsf.nb_representatives)
    found = 0
    start = time.time()
    for i in range(2**rsf.nb_representatives):
        rsf.set_SANF(tt.next_entry())
        if (resiliency==-1 or rsf.is_resilient_optimised(resiliency)) and rsf.is_algebraic_immune(algebraic_immunity):
            found += 1
            rsf.update_TT()
            print(rsf.TT)
           
    end = time.time()
    print("time elapsed: " + str(end - start) + " s")
    return found
    

   
def find_RSF(locality, resiliency, algebraic_immunity, max_degree_SANF = [], min_degree_SANF = [0]):
    """
    Exhaustive approach optimised for dahus.\n
    Resiliency and algebraic immunity must still be specified.\n
    The exhaustive search is made on the SANF and the representives of degree above (l+1)/2 are set to 0. \n
    The max_degree_SANF and min_degree_SANF parameters are usefull to parallelize the computation:\n
    -The max_degree_SANF parameter allows to constraint the SANF of the representives of maximal degree (l+1)/2.
    It must contain the exact number of such representatives.\n
    -The min_degree_SANF parameter allows to fix the minimal degree representatives of the SANF.
    It can contain any number of elements. By default, it is [0], meaning that the functions only searches RSF not containing +1.\n

    Results
    -----------------------
    The results are found in a file in the result directory and contain all SANF and ANF.\n
    The filename is of the form "rsf-<locality>-<resiliency>-<AI>.txt".\n
    This function creates a backup file in the backup directory, it is updated every 30 minutes.\n
    In case the function is aborted, a second call to the function will resume at the last backup.\n
    Because of this mechanism, it might be possible that a function appears twice in the result file.\n
    If the function finishes normally, an "End" tag ends to the result file.
    

    Parameters
    ----------
    locality : integer
        number of variables to consider.
    resiliency : integer
        resiliency to verify (should be optimal).
    algebraic_immunity : integer
        algebraic immunity to verify (should be optimal).
    max_degree_SANF : array of Booleans, optional
        SANF of maximal degree. Empty by default.
    min_degree_SANF : array of Booleans, optional
        SANF of small degrees, the array can contain any number of elements. The default is [0].


    """
    
    rsf = RSF(locality)
    
    #maximal degree of the dahu
    max_degree = int((locality+1)/2)
    print("Maximal degree: " + str(max_degree))
    
    #number of representatives of maximal degree
    nb_max_degree_representatives = 0
    while sum(rsf.representatives[nb_max_degree_representatives]) <= max_degree:
        nb_max_degree_representatives += 1
        
    #number of representatives of lower degrees
    nb_low_degree_representatives = 0
    while sum(rsf.representatives[nb_low_degree_representatives]) < max_degree:
        nb_low_degree_representatives += 1
    nb_max_degree_representatives = nb_max_degree_representatives - nb_low_degree_representatives
    
    print("Nb of representatives: " + str(rsf.nb_representatives))
    print("Nb of representatives of maximal degree: " + str(nb_max_degree_representatives))
    print("Nb of representatives of lower degrees: " + str(nb_low_degree_representatives))
    
    #fix the SANF of higher degrees to zero
    high_degree_SANF = [0]*(rsf.nb_representatives - nb_low_degree_representatives - nb_max_degree_representatives)
    
    #SANF of maximal degree
    #if not specified by user, exhaustive search.
    if len(max_degree_SANF) != nb_max_degree_representatives:
        max_degree_SANF = []
        nb_low_degree_representatives += nb_max_degree_representatives
        print("No representative of maximal degree is given or misformed: exhaustive search.")
    
    #exhaustive search on SANF of lower degrees, except those already fixed by user
    low_degree_SANF = Truth_table_entry(nb_low_degree_representatives - len(min_degree_SANF))
    
    #result and backupo files
    fichier_resultat = open("result/rsf-"+str(locality)+"-"+str(resiliency)+"-"+str(algebraic_immunity)+"-"+''.join([str(i) for i in max_degree_SANF])+"-"+''.join([str(i) for i in min_degree_SANF])+".txt","a")
    fichier_backup_name = "backup/rsf-"+str(locality)+"-"+str(resiliency)+"-"+str(algebraic_immunity)+"-"+''.join([str(i) for i in max_degree_SANF])+"-"+''.join([str(i) for i in min_degree_SANF])+".txt"
    
    #search and load backup
    try:
        fichier_backup = open(fichier_backup_name)
        tmp = locals()
        exec(fichier_backup.read(),globals(),tmp)
        fichier_backup.close()
        backup = tmp["backup"]
        if backup == 'End': #if the backup says the computation has ended
            print("Backup found.\nAll results are already computed!\nSee the result directory.")
            return 0
        low_degree_SANF.current = backup
        print("Backup found.")
        print("Backup: "+ str(low_degree_SANF.current))
        
    except IOError:
        print("No backup found.")
        
    
    fichier_backup = open(fichier_backup_name,"w")
    
    #compute how many SANF are still to be checked
    if 0 not in low_degree_SANF.current: #no backup
        stop = 2**(low_degree_SANF.locality)
    else: #backup
        stop = 2**(low_degree_SANF.locality) - bool_list_to_integer(low_degree_SANF.current)
    
    found = 0
    start = time.time()
    interval = time.time() -3600 #force backup at the very beginning
    
    for i in range(stop):
        
        #build SANF by concetenating every part, the exhaustive search is made on low_degree_SANF
        rsf.set_SANF(min_degree_SANF + low_degree_SANF.next_entry() + max_degree_SANF + high_degree_SANF)
        
        #verify resiliency and AI
        if rsf.is_resilient_optimised(resiliency) and rsf.is_algebraic_immune(algebraic_immunity):
            found += 1
            rsf.update_ANF_from_SANF()
            fichier_resultat.write(str(rsf.SANF) + "\n")
            fichier_resultat.write(str(rsf.ANF) + "\n")
            fichier_resultat.flush()
            
        #backup
        if time.time() - interval > 1800:
            interval = time.time()
            fichier_backup.write("backup=" + str(low_degree_SANF.current) + "\n")
            fichier_backup.flush()
            
    #fermeture
    end = time.time()
    fichier_resultat.write("End")
    fichier_resultat.flush()
    fichier_resultat.close()
    fichier_backup.write("backup='End'")
    fichier_backup.flush()
    fichier_backup.close()
    
    return found

def find_RSF_with_coverage(locality, resiliency, algebraic_immunity, max_degree_SANF, min_degree_SANF = [0], out=print):
    """
    Exhaustive approach optimised for dahus.\n
    Resiliency and algebraic immunity must still be specified.\n
    The exhaustive search is made on the SANF with some restrictions:\n
    -the representives of degree above (l+1)/2 are set to 0. \n
    -the representives of maximal degree (l+1)/2 are fixed with the max_degree_SANF parameter. \n
    -the representatives not covered by one of the set maximal degree representatives are set to 0.\n
    
    The max_degree_SANF and min_degree_SANF parameters are usefull to parallelize the computation:\n
    -The max_degree_SANF parameter allows to constraint the SANF of the representives of maximal degree (l+1)/2.
    It must contain the exact number of such representatives.\n
    -The min_degree_SANF parameter allows to fix the minimal degree representatives of the SANF.
    It can contain any number of elements. By default, it is [0], meaning that the functions only searches RSF not containing +1.\n

    Results
    -----------------------
    The results are found in a file in the result directory and contain all SANF and ANF.\n
    The filename is of the form "rsf-c-<locality>-<resiliency>-<AI>.txt".\n
    This function creates a backup file in the backup directory, it is updated every 30 minutes.\n
    In case the function is aborted, a second call to the function will resume at the last backup.\n
    Because of this mechanism, it might be possible that a function appears twice in the result file.\n
    If the function finishes normally, an "End" tag ends to the result file.
    
    Limitations
    -----------
    So far, this algorithm has been tested with a single maximal degree representative set to 1.

    Parameters
    ----------
    locality : integer
        number of variables to consider.
    resiliency : integer
        resiliency to verify (should be optimal).
    algebraic_immunity : integer
        algebraic immunity to verify (should be optimal).
    max_degree_SANF : array of Booleans
        SANF of maximal degree.
    min_degree_SANF : array of Booleans, optional
        SANF of small degrees, the array can contain any number of elements. The default is [0].
    out : function, optional
        display function, print by default.


    """
    
    out('Start process id: ', os.getpid())
    
    rsf = RSF(locality)
    
    #maximal degree of the dahu
    max_degree = int((locality+1)/2)
    out("Maximal degree: " + str(max_degree))
    
     #number of representatives of maximal degree
    nb_max_degree_representatives = 0
    while sum(rsf.representatives[nb_max_degree_representatives]) <= max_degree:
        nb_max_degree_representatives += 1
        
    #number of representatives of lower degrees
    nb_low_degree_representatives = 0
    while sum(rsf.representatives[nb_low_degree_representatives]) < max_degree:
        nb_low_degree_representatives += 1
    nb_max_degree_representatives = nb_max_degree_representatives - nb_low_degree_representatives
    
    out("Nb of representatives: " + str(rsf.nb_representatives))
    out("Nb of representatives of maximal degree: " + str(nb_max_degree_representatives))
    out("Nb of representatives of lower degrees: " + str(nb_low_degree_representatives))
    
    #fix the SANF of higher degrees to zero
    high_degree_SANF = [0]*(rsf.nb_representatives - nb_low_degree_representatives - nb_max_degree_representatives)
    
    #maximal degree SANF
    if len(max_degree_SANF) != nb_max_degree_representatives:
        out("Error: representatives of maximal degree must be specified!")
        return 0
    
    #indexes of maximal degree representatives set to 1
    max_degree_representative_indexes = []
    for i in range(len(max_degree_SANF)):
        if max_degree_SANF[i] == 1:
            max_degree_representative_indexes.append(i + nb_low_degree_representatives)
    
    
    #search representatives covered by maximal degree representatives set to 1
    covered_representatives = []
    for i in range(nb_low_degree_representatives):
        j = 0
        while j<len(max_degree_representative_indexes):
            if representative_cover(locality, max_degree_representative_indexes[j], i):
                covered_representatives.append(i)
                j = len(max_degree_representative_indexes)
            j+=1
    out("Nb of covered representatives: ", len(covered_representatives))
    
    #SANF of covered representatives
    covered_SANF = Truth_table_entry(len(covered_representatives) - len(min_degree_SANF))
    
    #initial SANF
    SANF = ([0]*nb_low_degree_representatives) + max_degree_SANF + high_degree_SANF
    for i in range(len(min_degree_SANF)):
        SANF[covered_representatives[i]] = min_degree_SANF[i]
    offset_index_coverage = len(min_degree_SANF)
    
    
    #result and backup files
    fichier_resultat = open("result/rsf-c-"+str(locality)+"-"+str(resiliency)+"-"+str(algebraic_immunity)+"-"+''.join([str(i) for i in max_degree_SANF])+"-"+''.join([str(i) for i in min_degree_SANF])+".txt","a")
    fichier_backup_name = "backup/rsf-c-"+str(locality)+"-"+str(resiliency)+"-"+str(algebraic_immunity)+"-"+''.join([str(i) for i in max_degree_SANF])+"-"+''.join([str(i) for i in min_degree_SANF])+".txt"
    
    #search and load an existing backup
    try:
        fichier_backup = open(fichier_backup_name)
        tmp = locals()
        exec(fichier_backup.read(),globals(),tmp)
        fichier_backup.close()
        backup = tmp["backup"]
        if backup == 'End':
            print("Backup found.\nAll results are already computed!\nSee the result directory.")
            return 0
        covered_SANF.current = backup
        out("Backup found.")
        out("Backup: "+ str(covered_SANF.current))
        
    except IOError:
        out("No Backup found.")
        
    
    fichier_backup = open(fichier_backup_name,"w")
    
    #compute how many SANF are still to be checked
    if 0 not in covered_SANF.current: #backup
        stop = 2**(covered_SANF.locality)
    else: #no backup
        stop = 2**(covered_SANF.locality) - bool_list_to_integer(covered_SANF.current)
    
    found = 0
    interval = time.time() -7200 #force a backup at the very beginning
    
    for i in range(stop):
        
        #combine SANF
        covered_SANF.next_entry()
        for i in range(covered_SANF.locality):
            SANF[covered_representatives[i + offset_index_coverage]] = covered_SANF.current[i]
        
        #set SANF
        rsf.set_SANF(SANF)
        
        #check resiliency and AI
        if rsf.is_resilient_optimised(resiliency) and rsf.is_algebraic_immune(algebraic_immunity):
            found += 1
            rsf.update_ANF_from_SANF()
            fichier_resultat.write(str(rsf.SANF) + "\n")
            fichier_resultat.write(str(rsf.ANF) + "\n")
            fichier_resultat.flush()
            
        #backup
        if time.time() - interval > 1800:
            interval = time.time()
            fichier_backup.write("backup=" + str(covered_SANF.current) + "\n")
            fichier_backup.flush()
            
    end = time.time()
    fichier_resultat.write("End")
    fichier_resultat.flush()
    fichier_resultat.close()
    fichier_backup.write("backup='End'")
    fichier_backup.flush()
    fichier_backup.close()
    
    return found