"""
Created on Sun Feb 28 11:47:53 2021

author: N. Abrate.

file: .py

description: Convert .txt with Marshak integrals (half-range, positive)
             produced with Mathematica (Byers formula) to .txt for TEST
"""
import numpy as np
from copy import deepcopy as cp

coeffs0 = np.loadtxt('MathematicaMarshak2000.txt')
NMax = int(np.sqrt(coeffs0.shape[0]))
coeffs0 = coeffs0[:, 2].reshape((NMax, NMax), order='C')
Mcoeffs = cp(coeffs0)
# coeffsN = np.loadtxt('Nicolo.txt')

PN = 2001
m = 1
M = PN
normcoeff = (2*(np.arange(0, PN+1))+1)/2

for row in range(0, (M+1)//2):
    for n in range(0, PN+1):
        if n % 2 == 0:
            pos = normcoeff[n]*coeffs0[m, n]
            Mcoeffs[row+(M+1)//2, n] = -pos
        else:
            if n == m:
                pos = 1/2
                Mcoeffs[row+(M+1)//2, n] = pos
            else:
                pos = 0
                Mcoeffs[row+(M+1)//2, n] = pos
        Mcoeffs[row, n] = pos
    m = m+2

np.savetxt('Marshak.txt', Mcoeffs[0:NMax//2, :], fmt='%.18e')
