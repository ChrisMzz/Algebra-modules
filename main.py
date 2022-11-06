from lib.source.equations import *
from lib.source.group import *
from lib.source.matrix import *

#print(SymmetricGroup(3).group)
#print(Equation(8).coeffs)
#print(SqMatrix(((1,2),(3,4))).caracteristic_p())

S5 = SymmetricGroup(5)

print("S5 Permutation, orbit, order, signature")
for element in S5:
    print(S5.to_string(element), S5.orbits(element), S5.order(element), S5.signature(element))


