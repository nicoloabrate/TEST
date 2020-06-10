"""
Created on Sat May 16 10:11:57 2020

author: N. Abrate.

file: .py

description:
"""
import time as t
from TEST.geometry import Slab
import TEST.NeutronTransportEquation as NTE
import TEST.eigenproblems.delta as delta
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm

nev = 40
M = 100
G = 1
N = 7
bc = 'Mark'
H = 10

matname = ['Modak']
xlayers = [0, H]
# define geometry and mesh
slab = Slab(M, xlayers, matname, [bc], G)
start = t.time()
PN = NTE.PN(slab, N, fmt='csc')
print('Elapsed time PN: %f' % (t.time()-start))

d = delta(slab, PN, nev=nev, algo='eigs', verbosity=True)
v = np.real(d.eigvect[0:slab.NT, 0])

fig = plt.figure()
plt.plot(slab.mesh, v/norm(v))
