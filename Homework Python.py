# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 15:19:01 2019

@author: m1xam
"""
from scipy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
import sys
import pylab
from scipy.integrate import quad
from scipy.optimize import fsolve
from math import sqrt

n=10
#Task 1
def approx_ln(x, n):

    a = (1 + x)/2
    g = sqrt(x)
    for n in range(n):
        a = (a+g)/2
        g = sqrt(a*g)
    return((x-1)/a)

print(approx_ln(100,10))
print(np.log(100))

# Task 2

A = []
B = []
C = []
def approx_ln(x, n):
    a = (1 + x)/2
    g = sqrt(x)
    for n in range(n):
        A.append((x-1)/a)
        B.append(np.log(x))
        C.append((x-1)/a-np.log(x))
        a = (a+g)/2
        g = sqrt(a*g)
   
approx_ln(1.41,20)
matplotlib.pyplot.figure(1)  
plot(A)
plot(B)

matplotlib.pyplot.figure(2)
plot(C)

### Task 3

C = []
def approx_ln(x, n):
    a = (1 + x)/2
    g = sqrt(x)
    for n in range(n):
        C.append(abs((x-1)/a-np.log(x)))
        a = (a+g)/2
        g = sqrt(a*g)

    
approx_ln(1.41,20)

matplotlib.pyplot.figure(3)
plot(C)
#
### Task 4
 
M = []
def fast_approx_ln(x, n):
    A = []
    a = (1 + x)/2
    g = sqrt(x)
    for k in range(n):
        A.append(a)
        a = (a+g)/2
        g = sqrt(a*g)
    M.append(A)

    for k in range(1, len(M[-1])):
        B=[]
        for e in range(len(M[-1]) - 1):
            d = (M[-1][e+1] - 4**(-k)*M[-1][e])/(1-4**(-k))
            B.append(d)
        M.append(B)
    return((x-1)/d)
    
print(fast_approx_ln(100,2))

# Task 5

M = []
def fast_approx_ln(x, n):
    A = []
    a = (1 + x)/2
    g = sqrt(x)
    for k in range(n):

        a = (a+g)/2
        g = sqrt(a*g)
        A.append(a)
    M.append(A)

    for k in range(1, len(M[-1])):
        B=[]
        for e in range(len(M[-1]) - 1):
            d = (M[-1][e+1] - 4**(-k)*M[-1][e])/(1-4**(-k))
            B.append(d)

        M.append(B)
        
    return((x-1)/d)
print(fast_approx_ln(100,10))
X_2 = []
X_3 = []
X_4 = []
X_5 = []

z = linspace(0.1,20,250)

for i in range(len(z)):
    x = z[i]
    X_2.append(abs(np.log(x)-fast_approx_ln(x,3)))
    X_3.append(abs(np.log(x)-fast_approx_ln(x,4)))
    X_4.append(abs(np.log(x)-fast_approx_ln(x,5)))
    X_5.append(abs(np.log(x)-fast_approx_ln(x,6)))

plt.figure(4)

plt.xlabel("x")
plt.ylabel("Error")
plt.title("Error behaviour of the accelerated Carlsson Method for the log")

plt.xlim(0,20)
plt.ylim(1.e-19,1.e-5)

plt.semilogy(z,X_2, label = "2 iterations")
plt.semilogy(z,X_3, label = "3 iterations")
plt.semilogy(z,X_4, label = "4 iterations")
plt.semilogy(z,X_5, label = "5 iterations")
legend(loc='upper left')



















  