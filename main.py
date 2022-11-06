#from lib.source.equations import *
from lib.source.group import *
from lib.source.matrix import *

#print(SymmetricGroup(3).group)
#print(Equation(8).coeffs)
#print(SqMatrix(((1,2),(3,4))).caracteristic_p())

A3 = AlternatingGroup(3)
Z2 = CyclicGroup(2)
K = KroneckerProduct(Z2, A3)

print("A3 x Z2, A3 element order, Z2 element order, A3 x Z2 element order")
for element in K.group:
    print(K.to_string(element), ",", K.group1.order(element[0]), ",", K.group2.order(element[1]), ",", K.order(element))



