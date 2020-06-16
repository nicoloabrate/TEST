"""
Created on Sat May 16 10:11:57 2020

author: N. Abrate.

file: .py

description:
"""
import time as t
from TEST.geometry import Slab
import TEST.NeutronTransportEquation as NTE
from TEST.eigenproblems.alpha import alphaprompt
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm


Modak_benchmark = [-2.53782e-2, -1.03353e-1, -2.38123e-1,
                   -4.47714e-1]

# FIXME: add S8 with 100 meshes and S16 with 50 meshes to have better comparison

nev = 5
M = 50
G = 1
N = 15
bc = 'Mark'
H = 10

matname = ['Dahl']
xlayers = [0, H]
# define geometry and mesh
slab = Slab(M, xlayers, matname, [bc], G)
start = t.time()
PN = NTE.PN(slab, N, fmt='csc', steady=False, prompt=True)
print('Elapsed time PN: %f' % (t.time()-start))

a = alphaprompt(slab, PN, nev=nev, algo='eigs', verbosity=True)
v = np.real(a.eigvect[0:slab.NT, 0])
