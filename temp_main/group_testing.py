from source.group import *
from source.equations import *
import numpy as np

np.random.seed(42)

#Exemples

# Créer des groupes :
#S = SymmetricGroup(4)
#A = AlternatingGroup(3)
#I = InverseCyclicGroup(10)
#C = CyclicGroup(4)



# Afficher les éléments :
#print(A.group)
#print(S.group)

# -------------------------------------------------------------------------------
# Pour exporter en .csv (et pouvoir ouvrir les résultats dans Excel)
# 
# print("colonne 1 ; colonne 2 ; etc.")     # première ligne
# print(données1_col1, ";", données1_col2, ";", etc)  # ligne 2
# print(données2_col1, ";", données2_col2, ";", etc)  # ligne 3
# ...

# Ouvrez un invité de commandes dans le dossier contenant ce fichier et tapez :
# python main.py > fichier.csv

# -------------------------------------------------------------------------------

"""
A3 = AlternatingGroup(3)
Z2 = CyclicGroup(2)
Z9x = UnitGroup(9)
K = KroneckerProduct(Z2, A3)

print("element ; order ; ; element ; order1 ; order 2 ; order")
for element in Z9x.group:
    print(Z9x.to_string(element), ";", Z9x.order(element))
print()
for element in K.group:
    print(K.to_string(element), ";", K.group1.order(element[0]), ";", K.group2.order(element[1]), ";", K.order(element))


EQ1 = Equation(8)
print(EQ1.to_string(0))
print(EQ1.to_string(1))
print(EQ1.solutions)

"""


S = SymmetricGroup(5)
print("Permutation ;", "Orbites ;", "Ordre ;", "Signature")
for perm in S.group:
    print(S.to_string(perm), " ; ", S.orbits(perm), " ; ", S.order(perm), " ; ", S.signature(perm))

