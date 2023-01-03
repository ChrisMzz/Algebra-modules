from algch_chris_mzz.group import*
from copy import deepcopy
import sympy as sp
import numpy as np


def sort(M, subindex_list):
    si = [i[0] + i[1] for i in subindex_list]
    dic = {}
    for i in range(len(si)):
        dic[si[i]] = i
    l = []
    for index in sorted(si):
        l.append(M[dic[index]])
    return l

class BinomialCoeff : 
    def __init__(self,n,p):
        self.n = n
        self.p = p
        self.card = (factorial(n)/(factorial(n-p)*factorial(p)))
        self.elements = []
        self.build()
        
    def __str__(self):
        temp = ""
        for element in self.elements:
            temp += str(element) + "\n"
        return temp
    
    def __iter__(self):
        for element in self.elements:
            yield element
    
    def build(self):
        for sigma in SymmetricGroup(self.n):
            pot = [int(el) for el in sigma]
            while len(pot) != self.p:
                pot.pop(-1)
            if sorted(pot) not in self.elements:
                self.elements.append(sorted(pot))
    
    
        

class Matrix :
    def __init__(self, *args) -> None:
        if type(args[0]) == tuple:
            line_list = args[0]
        else:
            line_list = args
        self.matrix = [list(line) for line in line_list]
        self.n = len(self.matrix[0])
        m = 0
        for line in self.matrix:
            if len(line) != self.n:
                self.matrix = [[]]
                self.n = 0
                print("Matrix lines do not align.")
                break
            m += 1
        self.m = m
        
    def __str__(self):
        lines = ""
        for line in self.matrix:
            lines += str(line) + "\n"
        return lines
        


class SqMatrix (Matrix):
    def __init__(self, *args) -> None:
        super().__init__(*args)
        if self.n != self.m:
            self.matrix = [[]]
            #print("Matrix does not have as many lines as columns.")
            self.n, self.m = 0, 0
        self.coeffs = self.characteristic_polynomial_coefficients()
        x = sp.Symbol("x")
        self.characteristic_polynomial = 0
        for i in range(self.n+1):
            self.characteristic_polynomial += self.coeffs[i]*x**(self.n-i)
        self.eigenvalues = self.get_eigenvalues()
    
    def __mul__(self, other):
        if type(other) == SqMatrix:
            return self.prod(other)
        else:
            A = self.matrix
            for i in range(self.n):
                for j in range(self.n):
                    A[i][j] *= other
            return SqMatrix(tuple(A))
                    
    def __add__(self, other):
        A = deepcopy(self.matrix)
        if type(other) != SqMatrix:
            temp_tuple = ()
            for i in range(self.n):
                line = []
                for j in range(self.n):
                    if i == j:
                        line.append(other)
                    else:
                        line.append(0)
                temp_tuple += (line,)
            B = SqMatrix(temp_tuple).matrix
        else:
            B = other.matrix
        for i in range(self.n):
            for j in range(self.n):
                A[i][j] += B[i][j]
        return SqMatrix(tuple(A))
    
    def __sub__(self, other):
        A = deepcopy(self.matrix)
        if type(other) != SqMatrix:
            temp_tuple = ()
            for i in range(self.n):
                line = []
                for j in range(self.n):
                    if i == j:
                        line.append(other)
                    else:
                        line.append(0)
                temp_tuple += (line,)
            B = SqMatrix(temp_tuple).matrix
        else:
            B = other.matrix
        for i in range(self.n):
            for j in range(self.n):
                A[i][j] -= B[i][j]
        return SqMatrix(tuple(A))
    
    def __pow__(self, k):
        if k == 0:
            temp_tuple = ()
            for i in range(self.n):
                line = []
                for j in range(self.n):
                    if i == j:
                        line.append(1)
                    else:
                        line.append(0)
                temp_tuple += (line,)
            return SqMatrix(temp_tuple)
        if k == 1:
            return self
        temp = self.prod(self)
        if k == 2:
            return temp        
        for i in range(k-2):
            temp = self.prod(temp)
        return temp
    
    def transpose(self):
        """Compute the transpose of the matrix.

        Returns:
            SqMatrix: The transpose of the matrix.
        """
        lines = {}
        n = self.n
        for i in range(n):
            lines[i] = []
        for i in range(n):
            for j in range(n):
                lines[j].append(self.matrix[i][j])
        return SqMatrix(tuple(lines.values()))
    
    def commutes_with(self, matrix):
        """Check whether the matrix commutes with another matrix.

        Args:
            matrix (SqMatrix): The other matrix.

        Returns:
            bool: Whether the matrices commute or not.
        """
        lr = self.prod(matrix)
        rl = matrix.prod(self)
        if lr.matrix == rl.matrix:
            return True
        else:
            return False
    
    
    def trace(self):
        """Compute the trace of the matrix.

        Returns:
            float: The trace.
        """
        trace = 0
        for i in range(self.n):
            for j in range(self.n):
                if i == j:
                    trace += self.matrix[i][j]
        return trace
    
    def det(self, show=False):
        """Compute the determinant of the matrix.

        Args:
            show (bool, optional): Whether to display the computation. Defaults to False.

        Returns:
            float: The determinant.
        """
        n = self.n
        temp_sum = 0
        temp_str = ""
        Sn = SymmetricGroup(n)
        for perm in Sn:
            temp_prod = 1
            temp_prod_str = ""
            for i in range(n):
                temp_prod_str += f"{self.matrix[i][int(perm[i])-1]}*"
                temp_prod *= self.matrix[i][int(perm[i])-1]
            temp_prod_str = temp_prod_str[:len(temp_prod_str)-1]
            signature = str(Sn.signature(perm))[:len(str(Sn.signature(perm)))-1] + " "
            if len(signature) == 1 and perm != Sn.group[0]:
                signature = "+ "
            elif perm == Sn.group[0]:
                signature = ""
            temp_str += f"{signature}" + temp_prod_str + " "
            temp_sum += Sn.signature(perm)*temp_prod
        if show:
            print(temp_str)
        return temp_sum    
    
    
    def prod(self, matrix):
        """Compute a matrix product.

        Args:
            matrix (SqMatrix): The other matrix.

        Returns:
            SqMatrix: _description_
        """
        n = self.n
        new_matrix = ()
        for i in range(n):
            temp_list = []
            for j in range(n):
                temp = 0
                for k in range(n):
                    temp += self.matrix[i][k]*matrix.matrix[k][j]
                temp_list.append(temp)
            new_matrix += (temp_list,)
        return SqMatrix(new_matrix)
    
    def vec_prod(self, vec):
        """Compute the application of a matrix to a vector.

        Args:
            vec (list): Vector (list of coordinates)

        Returns:
            list: The new vector.
        """
        n =self.n
        if len(vec) != n:
            return
        coords = []
        for i in range(n):
            tempsum = 0
            for j in range(n):
                tempsum += self.matrix[i][j]*vec[j]
            coords.append(tempsum)
        return tuple(coords)
    
    def characteristic_polynomial_coefficients(self):
        """Get characteristic polynomial coefficients in decreasing order.

        Returns:
            list: List of coefficients.
        """
        coeffs = [(-1)**(self.n), (-1)**(self.n-1)*self.trace()]
        for i in range(2,self.n): # fait varier le coefficient étudié, du coefficient du terme de plus haut degré au coefficient constant
            temp_c = 0
            for ij in BinomialCoeff(self.n,i):
                M = [[self.matrix[k-1][k-1]] for k in ij] # sera le mineur principal donné par rapport à J = ij
                index_list = [[(k-1,k-1)] for k in ij]
                for current_coeff in BinomialCoeff(i,2):
                    pos = [ij[k-1]-1 for k in current_coeff]
                    for line in range(len(M)): # nécessaire pour ordonner la matrice créée, et calculer le "bon" déterminant du mineur M
                        if (pos[0],pos[0]) in index_list[line]:
                            index_list[line].append((pos[0],pos[1]))
                            M[line].append(self.matrix[pos[0]][pos[1]])
                        if (pos[1],pos[1]) in index_list[line]:
                            index_list[line].append((pos[1],pos[0]))
                            M[line].append(self.matrix[pos[1]][pos[0]])
                for line in range(len(M)): # ordonne chaque ligne de M en fcontion de index_list qui dépend de ij
                    M[line] = sort(M[line],index_list[line])        
                temp_c += SqMatrix(tuple(M)).det() # calcule le déterminant du mineur M
            coeffs.append((-1)**(self.n-i)*temp_c) 
            # J'ai écrit n-i et non i car contrairement à la démonstration de cours j'étudie les termes en ordre décroissant
        coeffs.append(self.det())
        return coeffs
    
    
    def is_colinear_with_basis_element(self, vector, basis):
        """Check if a vector is colinear with an element from a basis.

        Args:
            vector (list): Vector (list of coordinates)
            basis (list): Basis (list of vectors)

        Returns:
            bool: Whether there exists a colinear vector in the basis or not.
        """
        n = len(vector)
        for basis_vector in basis:
            k = 0
            while basis_vector[k] == 0:
                k += 1
            if vector[k] == 0:
                pass
            else:
                (l, m) = (basis_vector[k], vector[k])
                scaled_basis_vector = [basis_vector[i]*m for i in range(n)]
                scaled_vector = [vector[i]*l for i in range(n)]
                if scaled_basis_vector == scaled_vector:
                    return True
        return False
    

    def get_eigenvalue_algmultiplicity(self, eigenvalue):
        """Compute the algebraic multiplicity of a specified eigenvalue.

        Args:
            eigenvalue (float): The specified eigenvalue.

        Returns:
            int: The eigenvalue's algebraic multiplicity.
        """
        algmultiplicity = 0
        x = sp.Symbol("x")
        n = self.n
        coeffs = self.coeffs
        polynomial = 0
        for i in range(n+1):
            polynomial += coeffs[i]*x**(n-i)
        while eigenvalue in sp.solvers.solve(polynomial):
            algmultiplicity += 1
            polynomial = 0
            temp = []
            for i in range(n+1-algmultiplicity):
                temp.append((n-i-(algmultiplicity-1))*coeffs[i])
            coeffs = temp
            for i in range(n+1-algmultiplicity):
                polynomial += coeffs[i]*x**(n-i-algmultiplicity)
        return algmultiplicity
            
        
    def get_eigenvalues(self):
        """Compute eigenvalues associated to the matrix.

        Returns:
            dict: Dictionary of the form {eigenvalue : algebraic multiplicity}.
        """
        x = sp.Symbol("x")
        n = self.n
        coeffs = self.coeffs
        polynomial = 0
        for i in range(n+1):
            polynomial += coeffs[i]*x**(n-i)
        eigenvalues = sp.solvers.solve(polynomial, x)
        values_with_algmultiplicity = {}
        for eigenvalue in eigenvalues:
            values_with_algmultiplicity[eigenvalue] = self.get_eigenvalue_algmultiplicity(eigenvalue)
        return values_with_algmultiplicity

        
        
    def get_eigenvectors(self, eigenvalue):
        if eigenvalue not in self.eigenvalues.keys():
            return
        A = deepcopy(self)
        n = self.n
        x = []
        for i in range(1,n+1):
            x.append(sp.Symbol(f"x{i}"))
        return sp.solvers.solve((A-eigenvalue).vec_prod(x), x)
    
    def get_generalised_eigenvectors(self, eigenvalue):
        """Get all generalised eigenvectors associated to a specified eigenvalue

        Args:
            eigenvalue (float): The eigenvector

        Returns:
            list: List of generalised eigenvectors.
        """
        general_vectors = [self.get_eigenvectors(eigenvalue)]
        alg, p = self.eigenvalues[eigenvalue] - 1, 1
        A = deepcopy(self)
        n = self.n
        x = []
        for i in range(1,n+1):
            x.append(sp.Symbol(f"x{i}"))
        while alg > 0:
            p += 1
            general_vectors.append(sp.solvers.solve(((A-eigenvalue)**p).vec_prod(x), x))
            alg -= 1
        return general_vectors

    def make_vector_from_solution(self, solution, symbols, bias=0):
        """Transform a sympy dictionary solution into a vector.

        Args:
            solution (dict): The solution. Can also be an empty list.
            symbols (list): List of sympy symbols
            bias (int, optional): Shifts parameter selection (only used if previous selection is colinear with another vector)
        Returns:
            list: The vector created.
        """
        vector = []
        parameters = {}
        if bias >= self.n or solution == {}:
            return
        if solution == [] or solution == {}: # if ker(A-eigenvalue) = R^n
            pot = [[int(i==j) for i in range(self.n)] for j in range(self.n)]
            for symbol in symbols:
                f = sp.lambdify([symbols], symbol, 'numpy')
                vector.append(f(pot[bias]))
            return vector
        i = 0
        for symbol in symbols:
            if symbol not in solution.keys():
                parameters[symbol] = int(i==bias)
                i += 1
        #print(f"parameters for {solution} : {parameters}")
        for symbol in symbols:
            if symbol in solution.keys():
                f = sp.lambdify([parameters.keys()], solution[symbol], 'numpy')
                vector.append(f(list(parameters.values())))
            else:
                f = sp.lambdify([parameters.keys()], symbol, 'numpy')
                space_size = len(parameters.keys())
                pot = [[int(i==j) for i in range(space_size)] for j in range(space_size)]
                vector.append(f(pot[bias]))
        return vector
      
    def get_all_generalised_eigenvectors(self):
        """Get every generalised eigenvector, for each eigenvalue.

        Returns:
            list: List of eigenvectors.
        """
        space_basis = [] # list to be filled with generalised eigenvectors
        n = self.n
        x = []
        nul_vector = [0 for _ in range(n)]
        for i in range(1,n+1):
            x.append(sp.Symbol(f"x{i}"))
        for eigenvalue in self.eigenvalues.keys():
            if self.eigenvalues[eigenvalue] == 1:
                solution = self.get_eigenvectors(eigenvalue)
                print(f"Solution for {eigenvalue} : {solution}")
                space_basis.append(self.make_vector_from_solution(solution, x))
            else:
                for solution in self.get_generalised_eigenvectors(eigenvalue):
                    print(f"Solution for {eigenvalue} : {solution}")
                    bias=0
                    vector = self.make_vector_from_solution(solution, x, bias)
                    while vector in space_basis and vector != nul_vector and self.is_colinear_with_basis_element(vector, space_basis):
                        bias +=1
                        vector = self.make_vector_from_solution(solution, x, bias)
                    space_basis.append(self.make_vector_from_solution(solution, x, bias))
        return space_basis

    def triangular_with_general_eigenspaces(self, display=True, round=-1):
        """Get triangular matrix and change-of-basis made with generalised eigenvectors.

        Args:
            display (bool, optional): Whether to display results. Defaults to True.
            round (int, optional): Rounds all matrix elements to `round` decimal places.

        Returns:
            (SqMatrix, SqMatrix): (T, P), the triangular matrix and change-of-basis matrix.
        """
        n = self.n
        basis = self.get_all_generalised_eigenvectors()
        P = SqMatrix(tuple(basis)).transpose()
        inverse = np.linalg.inv(np.array(P.matrix))
        iP = []
        for i in range(n):
            temp = []
            for j in range(n):
                    temp.append(inverse[i][j])
            iP.append(temp)
        iP = SqMatrix(tuple(iP))
        T = iP*self*P
        if round != -1:
            for i in range(n):
                for j in range(n):
                    if round == 0:
                        T.matrix[i][j] = int(np.round(T.matrix[i][j],round))
                    else:
                        T.matrix[i][j] = np.round(T.matrix[i][j],round)
        if display:
            print(f"T =\n{T}\nP =\n{P}\n A = PTP^(-1)")
        return T, P
        
        
            
