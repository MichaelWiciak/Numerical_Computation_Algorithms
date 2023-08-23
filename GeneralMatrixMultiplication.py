#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 08:54:35 2022

@author: michaelwiciak
"""
import random
import matplotlib.pyplot as plt
import time


def generalMatrixMultiplication(X, Y, result):
    
    # Must X Y are the matrixes to be multiplied together
    # Note: the resulting array need to be provided .
    # iterate through rows of X
    for i in range(len(X)):
       # iterate through columns of Y
       for j in range(len(Y[0])):
           # iterate through rows of Y
           for k in range(len(Y)):
               result[i][j] += X[i][k] * Y[k][j]
    
    return result

def differenceD(D,Array):
    n = len(Array)
    j=0
    for i in range(n-1):
        while (j<n-1 and Array[i]-Array[j]>D):
            j+=1
        if Array[i]-Array[j]==D:
            return 1
    return 0
 
def bubbleSort(A):
    #A is array of n size
    n = len(A)
    for i in range(n-2):
        for j in range (n-2-i):
            if A[j]>A[j+1]:
                tmp = A[j]
                A[j] =A[j+1]
                A[j+1] = tmp
    
    return A
           
def graphTest():
    Array = []
    x=[]
    y=[]
    for j in range(1, 1000000, 100000):
        for i in range(1, j, 1):
            Array.append(random.randint(1,10000))
            Array.sort()
            t0 = time.time()
            differenceD(random.randint(1, 100), Array)
            processTime = time.time() - t0
            y.append(processTime)
            x.append(len(Array))

    plt.plot(x, y) 
    plt.xlabel('Array size')
    # naming the y axis
    plt.ylabel('Time Taken')
      
    # giving a title to my graph
    plt.title('DifferenceD')
      
    # function to show the plot
    plt.show()           
    
graphTest()
            
        