from lib.source.equations import *
from lib.source.group import *
from lib.source.matrix import *
import time
#print(SymmetricGroup(3).group)
#print(Equation(8).coeffs)
#print(SqMatrix(((1,2),(3,4))).caracteristic_p())


A = SqMatrix(((2,-1,-1),(2,1,-2),(3,-1,-2)))
#A = SqMatrix(((7,4,0,0),(-12,-7,0,0),(20,11,-6,-12),(-12,-6,6,11)))


print(f"A = \n{A}")
print(A.characteristic_polynomial)
print(A.eigenvalues)
t = time.time()
A.triangular_with_general_eigenspaces(True,0)
print(f"Done in {time.time() - t} seconds.")
