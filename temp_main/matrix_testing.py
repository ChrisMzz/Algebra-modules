import numpy as np
from source.matrix import *




n = 5
tuple_a = tuple([tuple([np.random.randint(0,6) for i in range(n)]) for j in range(n)])
A = SqMatrix(tuple_a)


print(f"For wolfram :\n{tuple_a}\n")
print(A)
print(A.caracteristic_p())




