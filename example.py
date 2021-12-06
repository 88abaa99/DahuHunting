# -*- coding: utf-8 -*-
"""
@author: 
    
This test file replays some results exhibited in Section 7 of Eurocrypt Submission 63
"""

from BF import BF
from RSF import RSF

""" Verify the example of Dahu_9 given in Section 7.4 """
print("*********** Verifying results of Section 7.4 using the BF class ***********")
locality = 9
f = BF(locality) # Creation of the object Boolean Function of locality 9
truthTable_hexa = 0x69c3e14be916349ef8c3163c1e25c3e9aa95a55a167c4fb007b85d66e15eb88399999666c9666399073c6eb434fa9e41556a9a9536a66c69f84762a9cb81915f
tt = [] # Conversion to Boolean list
while truthTable_hexa > 0:
    tt = [truthTable_hexa&1] + tt
    truthTable_hexa >>= 1
tt = [0]*(2**locality - len(tt)) + tt
print("Truth table: ", tt)
f.set_TT(tt) # Set the truth table of the Boolean Function
f.update_ANF() # Update the ANF representation
print("ANF: ", f.string_ANF())
f.update_WS() # Update the Walsh Spectrum (necessary for the resiliency)
print("Algebraic immunity 5: ", f.is_algebraic_immune(5))
print("Resiliency 3: ", f.is_algebraic_immune(3))

print("\n*********** Verifying results of Section 7.4 using the RSF class ***********")
locality = 11
sanf = [0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
rsf = RSF(locality) #Creation of the objest Rotational Symmetric Function of locality 11
rsf.set_SANF(sanf)
rsf.update_ANF_from_SANF()
print("ANF: ", rsf.ANF)
rsf.update_STT()
rsf.update_SWS()
print("Algebraic immunity 6: ", rsf.is_algebraic_immune(6))
print("Resiliency 4: ", rsf.is_resilient(4))



from find_BF import find_BF_naive
from find_RSF import find_RSF

print("\n*********** Counting the number of elements in Dahu_4 ***********")
print("Truth tables:")
n = find_BF_naive(4,1,2)
print(n, "dahus found")
print("\n*********** Counting the number of RSF in Dahu_7 ***********")
n = find_RSF(7,2,4, min_degree_SANF=[])
print(n, "dahus found")