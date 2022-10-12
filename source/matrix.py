from group import*
from copy import deepcopy


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
        self.matrix = [line for line in line_list]
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
    
    def commutes_with(self, matrix):
        lr = self.prod(matrix)
        rl = matrix.prod(self)
        if lr.matrix == rl.matrix:
            return True
        else:
            return False
    
    
    def trace(self):
        trace = 0
        for i in range(self.n):
            for j in range(self.n):
                if i == j:
                    trace += self.matrix[i][j]
        return trace
    
    def det(self, show=False):
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
    
    def caracteristic_p(self):
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

              
                
