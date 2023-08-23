#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 11:05:42 2022

@author: michaelwiciak
"""

import numpy as np


def lower_triangular_solve(A, b):
    """
    Solve the system  A x = b  where A is assumed to be lower triangular,
    i.e. A(i,j) = 0 for j > i, and the diagonal is assumed to be nonzero,
    i.e. A(i,i) != 0.
    
    The code checks that A is lower triangular and converts A and b to
    double precision before computing.

    ARGUMENTS:  A   lower triangular n x n array
                b   right hand side column n-vector

    RETURNS:    x   column n-vector solution
    """



def upper_triangular_solve(A, b):
    """
    Solve the system  A x = b  where A is assumed to be upper triangular,
    i.e. A(i,j) = 0 for j < i, and the diagonal is assumed to be nonzero,
    i.e. A(i,i) != 0.
    
    The code checks that A is upper triangular and converts A and b to
    double precision before computing.

    ARGUMENTS:  A   upper triangular n x n array
                b   right hand side column n-vector

    RETURNS:    x   column n-vector solution
    """

    A = A.astype(np.float64) 
    b = b.astype(np.float64) #we need to store both arrays as floats.
    n = len(b) #we assign n to be the size of b (alternatively n,n=A.shape)
    
    for i in range(n): #i is both the row number and the index of the diagonal element
        for j in range (0,i): #lower was (i+1,n), for upper we cycle through first i values in the row (up to the diagonal)
            if np.abs(A[i,j]>1.e-12): 
                print('A is not upper triangular')
                return #if any of the values above the diagonal are non-zero, give up
            else:
                A[j,i] = A[j,i]/A[j,j] #refactor the coefficents of the matrix so that the leading diagonals are 1 and other elements are scaled
               
    for i in range(n): #we cycle through the diagonal elements
        if abs(A[i,i]<1.e-12): #checks if the leading diagonal is zero, give up if it is
            print(f"A({i},{i}) is zero")
            return
        else:
            b[i] = b[i]/A[i,i] #we need to scale the b vector to match our vector
            A[i,i] = 1. #we now force the diagonal elements equal to 1
           
    x = np.zeros([n,1]) #we need to initialise the x vector, this could be any array(n,1) as we will be overwriting these values anyway
   
    x[n-1] = b[n-1]
    for i in range (2,n+1): #we want to cycle through the n-i non-zero elements
        x[n-i] = b[n-i] #This is the term in front of the sum
        for j in range(n-i+1, n): #our j index is from the diagonal to the end of the row
            x[n-i] = x[n-i] - A[n-i,j]*x[j] #this is the sum term, we take off A[i,j]x[j] each cycle
    return x
               
#these are the values from question 3b WS03, after we have converted to upper triangular   
#A = np.array([[4,3,2,1],[0,5/4,3/2,7/4],[0,0,11/5,-3/5],[0,0,0,40/11]])
#b = np.array([[10],[9/2],[8/5],[40/11]])
#upper_triangular_solve(A, b) 

    
    
    
def gauss_elimination(A, b, verbose=False):
    """
    Reduce the system  A x = b  to upper triangular form, assuming that
    the diagonal is nonzero, i.e. A(i,i) != 0.
    
    Before computing A and b are converted to double precision.

    ARGUMENTS:  A   n x n matrix
                b   right hand side column n-vector

                verbose  (optional) if true print elimination steps

    RETURNS:    A   upper triangular n x n matrix
                b   modified column n-vector
    """
    A = A.astype(np.float64) 
    b = b.astype(np.float64) #we need to store both arrays as floats.
    n = len(b) #we assign n to be the size of b (alternatively n,n=A.shape)        
  
    for i in range(n): #cycling through all the diagonals
        if np.abs(A[i,i]<1.e-15): #we check the diagonals are non-zero
            print('A is singular!')
            return #if any are zero, give up

        # Gaussian elimination
        b[i] = b[i]/A[i,i] #divide each of our b elements by the corresponding diagonal
        A[i,:] = A[i,:]/A[i,i] #divide each of our rows by the diagonal value, this makes our diagonals equal to 1
        
        #now we make our matrix upper triangular
        for j in range(i+1,n): #we act on one row at a time
            temp=A[j,i] #this temporary variable is the factor we use in the row operations (this would be the ratio between values if our diag wasn't 1)
            b[j] = b[j]-b[i]*temp #current row's b, updated by taking away a factor of the higher row b
            for k in range(i,n): #across all values of the current row, perform the row operation 
                A[j,k] = A[j,k]-A[i,k]*temp #this is us performing the row operation on all the values in the row

    
    #create a new array to store the results
    x = np.zeros(n)  # or    x=b
    
    #this section of code should look exactly like the upper_triangular_solve code
    x[n-1] = b[n-1] 
    for i in range(2,n+1): 
        x[n-i] = b[n-i]
        for j in range(n-i+1, n):
            x[n-i] = x[n-i] - A[n-i,j]*x[j]
        
    return x

#these are the values from question 3b WS03    
#A = np.array([[4,3,2,1],[1,2,2,2],[1,1,3,0],[2,1,2,3]])
#b = np.array([[10],[7],[5],[8]])
#gaussian_elimination(A, b)  

    
    
def generate_test_set(n):
    """
    Generate a random solvable system of linear equations.

    ARGUMENTS:  n   size of the problem
                
    RETURNS:    A   n x n matrix for system of linear equations
                x   n vector of soltuions
                b   n vector for right hand side
    """
    
    B = np.random.rand(n, n)
    eps = 0.1
    
    A = eps * np.eye(n) + B * B.T
    x = np.ones((n, 1))
    b = np.matmul(A, x)
    
    return A, x, b