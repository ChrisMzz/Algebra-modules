from lib.source.group import *
import numpy as np

class Divisors :
    def __init__(self, n):
        self.group = [i for i in range(1, n) if gcd(i, n) != 1]

    def to_string(self, el):
        return "[" + str(el) + "]"


class Equation :
    
    def __init__(self, n):
        self.n = n
        self.coeffs = []
        self.roots = ()
        self.solutions = []
        self.compute_rnd_coeffs()
        
    def to_string(self, i):
        coeffs = self.coeffs
        if self.roots != ():
            a = self.roots[0]
            b = self.roots[1]
        else:
            return
        
        if coeffs[1] < 0:
            middle = f"- {abs(coeffs[1])}"
        else:
            middle = f"+ {abs(coeffs[1])}"
        last = ""
        if coeffs[2] < 0:
            last = f" - {abs(coeffs[2])}"
        elif coeffs[2] > 0:
            last = f" + {abs(coeffs[2])}"
        
        if a < 0:
            fact1 = f"(x + {abs(a)})"
        elif a > 0:
            fact1 = f"(x - {abs(a)})"
        else:
            fact1 = "x"
        if b < 0:
            fact2 = f"(x + {abs(b)})"
        elif b > 0:
            fact2 = f"(x - {abs(b)})"
        else:
            fact2 = "x"
        
        return (f"x^2 {middle}x{last} = 0  (mod {self.n})",f"{fact1}{fact2} = 0  (mod {self.n})")[i]
    
    def compute_rnd_coeffs(self):
        n = self.n
        ablist = [(a,b) for a in range(1,n) for b in range(1,n)]
        temp_solutions = []
        while temp_solutions == [] and ablist != []:
            (a,b) = ablist[np.random.randint(0,len(ablist))]
            #print((a,b))
            for d in Divisors(n).group:
                if d*(d+a-b) % n == 0:
                    temp_solutions.append(a)
                    temp_solutions.append(b)
                    temp_solutions.append(d+a)
                    
                elif d*(d+b-a) % n == 0:
                    temp_solutions.append(a)
                    temp_solutions.append(b)
                    temp_solutions.append(d+b)
            ablist.remove((a,b))
        if temp_solutions == []:
            print("No valid (a,b) tuple found.")
            return False
        i,j = (np.random.randint(-1,2), np.random.randint(-1,2))
        a += n*i
        b += n*j
        solutions = []
        for element in temp_solutions:
            element = element % n
            if element != 0 and element not in solutions:
                solutions.append(element)
        if a*b % n == 0:
            solutions.append(0)
        self.coeffs = [1, -(a+b), a*b]
        self.roots = (a,b)
        self.solutions = solutions
        return True



