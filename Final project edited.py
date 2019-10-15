from  scipy import *
from  matplotlib.pyplot import *
from numpy import *
from scipy.sparse import *
import unittest
import timeit
import timeit as ti
from numpy import dot, array #for testing 

# Task 1 and 6

# Here we create our class, and we have both extended by an extra parameter tol. Furthermore, we've extended the repr-
# method, which allows us to represent the original matrix in CSC form aswell. However, the CSR to CSC conversion takes
# place in Task 4.


class SparseMatrix:

    def __init__(self, matrix, tol=0):
        self.matrix=np.asarray(matrix)
        self.intern_represent = 'CSR'
        self.number_of_nonzeros = self.get_nr_non_zeros()
        self.tol = tol
        self.CSR_representation = [self.A(), self.IA(), self.JA()]
        self.matrix_width = len(self.matrix[0])
        self.CSC_representation = [self.convert_to_CSC()]
        
    def __repr__(self):
        representation = ('\nCSR_representation:\n' + 'A = '+ str(self.CSR_representation[0]) +
                          '\nIA = ' + str(self.CSR_representation[1]) +
                          '\nJA = '+ str(self.CSR_representation[2]) +
                          '\n')
        try:
            return representation + ('\nCSC_representation:\n' + str(self.CSC_representation[0]))
        except AttributeError:
            return representation


    def A(self):
        array = []
        
        for row in self.matrix:
            for number in row:
                if abs(number) > self.tol:
                    array.append(number)
        return array
        
    def IA(self):
        array = [0]
        
        non_zeros = 0
        for row in self.matrix:
            for number in row:
                if abs(number) > self.tol:
                    non_zeros += 1
            array.append(non_zeros)
        return array
    
    def JA(self):
        array = []
    
        for row in self.matrix:
            index = 0
            for number in row:
                if abs(number) > self.tol:
                    array.append(index)
                index += 1
        return array
    
#TASK 2

#Since our matrix has two dimensions, the np.nonzero-command gives a tuple
#containing two arrays. Either one of these will get the correct number of non-
#zeros. In other words, wether we choose [0] ,as in the code, or [1] is arbitrary
    
    def get_nr_non_zeros(self):
        return len(np.nonzero(self.matrix)[0])
    
# TASK 3

# We construct a method that can change any given element in the matrix.
# After having changed a given element we update our number of non-zeros.

    def change_element (self,i,j,new):
        self.matrix [i,j]=new
        self.number_of_nonzeros = self.get_nr_non_zeros()
        self.CSR_representation = [self.A(), self.IA(), self.JA()]
        return self.matrix
        
#TASK 4

#Create a method changing the internal representationen from CSR to CSC.

# An important observation to make is that when a matrix is converted into CSR
# most, but not all information is kept. The information that is lost, is how
# many (if any) columns with exclusively zeros there is to the right of the last column
# contating non-zeros. In order to do a conversion from CSR to CSC we also
# this information: len(self.matrix[0])
        
# both A and JA has one element per non-zero. We first check the corresponding
# indexes for the numbers in increasing order in JA. We then append the values
# of the same indexes from A into CSC_A. 
        
# To generate CSC_IA we observe that it will be the accumulated number of 
#non-zeros, counting column-wise. Combining this with the information that the 
#0's in JA is the value of the element of index 1 in CSC_IA, the number of 0's 
#and 1's in JA is the value of the element of index 2 in CSC_IA, we have all
# information that is needed. 
#Note that here the information about the width of the matrix is crucial. For 
#every column after the last non-zero, we should get the same number as for the
# last non-zero column.
        
# In order to generate CSC_JA we need to sort the non-zeros in the following
# order: Up to down, left to right. We begin by finding the 0's in JA. Since
# we care which place, and not which index, this 0 has, we add 1. Since the 
#numbers in JA are ordered from left to right, up and down, we can combine this
# information with IA in the following manner: We see how many rows we have to
# "jump down" in order for the first 0 in JA to be included. The first element
# in IA is always an arbitrary 0. We therefore subtract 1, and append however
# many rows it takes to reach this 0 into CSC_JA. Then we check for more 0's. 
#After that we search for 1's until and repeat until we've gone through the 
#entirety of JA. Note that by doing this, we sort the 0's (which are in the 
#0'th column of the matrix) upwards and down, then move right to the 1's in the
# next column etc. This is the way CSC_JA should be ordered. We repeat this 
# until we've appended everything into CSC_JA. 
    
    def convert_to_CSC(self):
        
        CSC_A = []
        for n in range(len(self.matrix[0])):
            for i in range(len(self.CSR_representation[2])):
                if self.CSR_representation[2][i] == n:
                    CSC_A.append(self.CSR_representation[0][i])

        CSC_IA = [0]
        for n in range(len(self.matrix[0])):
            CSC_IA.append(self.CSR_representation[2].count(n) + CSC_IA[-1])
        
        CSC_JA = []
        for n in range(len(self.matrix[0])):
            for i in range(len(self.CSR_representation[2])):
                if self.CSR_representation[2][i] == n:

                    num = i + 1
                    
                    for k in self.CSR_representation[1]:
                        if k >= num:
                            CSC_JA.append(self.CSR_representation[1].index(k) - 1)
                            break
                        
        self.intern_represent = 'CSC'
        
        self.CSC_representation = [CSC_A, CSC_IA, CSC_JA]
        return self.CSC_representation
# TASK 5

# Write a method checking if two such matrices are (exactly) equal.

# Both CSC and CSR start in the upper left corner, but CSC moves up to down,
# then left to right. CSR moves left to right, then up and down. To check
# if our generated CSC and CSR represent the same original matrix we do the
# following: If CSR for a given matrix is identical to CSC for the same matrix,
# but transposed, and vice versa we conclude that our CSR and CSC represent
# the same matrix.
        
    def CSR_equal_to_CSC(self, other):
        if  self.CSR_representation == other.CSC_representation and self.CSC_representation == other.CSR_representation:
            return (True)
        else: 
            return (False)   

# task 7 

# Addition of two sparse matrices. 

# Here we only had partial success. We were able to generate both IA and JA correctly, but weren't able to do so for A.
# We begin by checking if the matrices are of the same dimensions. If not, a message is returned saying that they aren't
# of the same dimensions.

# We divide the task into two cases: For a given row, either the first matrix has the most non-zero elements, or they have
# exactly the same non-zeros, or the other matrix has the most non-zero elements. We treat the first two in "if x>=y"
# and the third case in "if y > x".

# Recall that everything is done row-wise. To create JA_new we do the following:

# We append all elements from the first row of the matrix with the biggest number of non-zeros. We then check if there
# exists a non-zero in the other matrix has a different index(this means different value in JA). If it does, we append
# it aswell.

# We repeat this process for every row in the matrices. Note that appending in this order isn't how we want the end 
# result to be. 
# We want the ordering of JA to be from left-to-right since we're working with CSR. We therefore sort JA_new elements, 
# for each corresponding row in the original matrix, at the end of this task.

# In order to generate IA_new we create a temporary list temp. Very similarily to how we generated JA_new, we append
# all of the unique indexes for the non-zeros for each corresponding row in the original matrices. But, we append
# this into temp, and we're not concerned about the values of the elements, only how many elements there is. This is
# the exact information that IA holds in CSR form, how many non-zero elements there are per row. However, IA should
# start from the first row and work downwards. The values in IA should be the number of non-zeros in the current row
# + the sum of non-zeros above. 

# Since IA always should start with a 0 we define IA_new as a list with a zero to begin with. We then append a new
# element, which's value is the length of temp (numbers of non-zeros in the current row) + the latest element in IA_new
# (since the non-zeros should accumulate in IA).

# Onto A: We could solve on paper, but couldn't translate it into code. We tried a few different approaches and didn't
# quite manage to get the correct answer. However, we have still decided to show what we managed to do.

# Here's a brief overview of what we want to do. We check if both matrices has non-zeros in a given row. If they do,
# we check if they overlap. If they do, we should add the value of them. If not, the new matrix should have the non-
# zero value in that place. 

# The first approach we tried was to let A_new = self.CSR_representation[0] (which is the first matrix's A), and insert
# the unique non-zero elements from the other matrix (and add the overlapping non-zeros).
# The problems that we faced was mainly due to the length of the list being changed throughout the process, so accessing
# the elements we wanted was tough.

# The approach that we got the furthest with was by the following approach: We begin by checking which corresponding
# matrix has the biggest number of non-zeros. We append these into A_new. If the other matrix has any non-zeros in a
# different index we append these aswell. We repeat this until we've gone through every row.

# Here we have two crucial problems that we were unable to solve: If for instance we have the the first row of our first
# matrix is [0,1,5,7] and the first row of the other matrix is [2,0,0,0] we would append 1,5 and 7 first, then 2, because
# the first matrix's row has 3 non-zero's, and the second only has one. We would get [1,5,7,2], but we should get
# [2,1,5,7].
#Secondly, if there is an overlap we should sum the non-zeros together. We didn't manage to do so, we added the wrong
# elements together, because we didn't know how to access exactly the elements we wanted.

# As the code is right now, we generate one non-zero into A_new for every place there should be a non-zero. But the
# ordering isn't right row-wise. And we don't have code for adding elements together, so if there is an overlap, only
# the value of one of the elements is shown.




    def addition(self,other):
        
        if len(self.matrix[0]) == len(other.matrix[0]) and len(self.CSR_representation[1]) == len(other.CSR_representation[1]):
            pass
        else: 
            return 'The two CSRs must be of two matrices with the same dimensions'
        
        A_new = []
        IA_new = [0]
        JA_new = []
        for i in range(len(self.CSR_representation[1])-1):
            
            x = self.CSR_representation[1][i+1] - self.CSR_representation[1][i]
            y = other.CSR_representation[1][i+1] - other.CSR_representation[1][i]
            
            if x >= y:
            
                temp = []
                k = self.CSR_representation[1][i]
                for a in range(x):
                    temp.append(self.CSR_representation[2][k])
                    JA_new.append(self.CSR_representation[2][k])
                    A_new.append(self.CSR_representation[0][k])

                    k = k + 1
                    
               

                K = other.CSR_representation[1][i]
                for A in range(y):
                    
                    if other.CSR_representation[2][K+A] not in temp:
                        
                        temp.append(other.CSR_representation[2][K+A])
                        JA_new.append(other.CSR_representation[2][K+A])
                        A_new.append(other.CSR_representation[0][K+A])
                        
                                          
                IA_new.append(len(temp)+IA_new[-1])
                
            if y > x:
                temp = []
                k = other.CSR_representation[1][i]
                for a in range(y):
                    temp.append(other.CSR_representation[2][k])
                    JA_new.append(other.CSR_representation[2][k])
                    A_new.append(other.CSR_representation[0][k])
                    k = k + 1
    
                    
                K = self.CSR_representation[1][i]
                for A in range(x):
                    
                    
                    if (self.CSR_representation[2][K+A]) not in temp:
                        temp[0:len(temp)] = sorted(temp[0:len(temp)])
                        temp.append(self.CSR_representation[2][K+A])
                        JA_new.append(self.CSR_representation[2][K+A])
                        A_new.append(self.CSR_representation[0][K+A])
                        
                IA_new.append(len(temp)+IA_new[-1])             
            
            
            JA_new[IA_new[i]:IA_new[i+1]] = sorted(JA_new[IA_new[i]:IA_new[i+1]]) # Here we sort the order in JA_new
            
        return [A_new, IA_new, JA_new]

a = SparseMatrix([[0,1,0,0,0],
                  [0,0,0,0,0],
                  [0,0,2,0,3],
                  [4,5,0,0,6],
                  [0,0,0,0,0]])

b = SparseMatrix([[0,0,0,0,8],
                  [0,0,0,9,15],
                  [0,10,0,0,0],
                  [0,0,11,0,0],
                  [12,0,13,14,0]])