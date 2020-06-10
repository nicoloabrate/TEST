"""
Created on Sat May 16 10:11:57 2020

author: N. Abrate.

file: .py

description:
"""
import time as t
from test.geometry import Slab
import test.NeutronTransportEquation as NTE
import test.eigenproblems.gamma as gamma
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm

nev = 10
M = 100
G = 1
N = 63
bc = 'Mark'
H = 8

matname = ['Modak']
xlayers = [0, H]
# define geometry and mesh
slab = Slab(M, xlayers, matname, [bc], G)
start = t.time()
PN = NTE.PN(slab, N, fmt='csc')
print('Elapsed time PN: %f' % (t.time()-start))
g1 = gamma(slab, PN, nev=nev, verbosity=True)
v1 = np.real(g1.eigvect[0:slab.NT, 0])

g2 = gamma(slab, PN, nev=nev, algo='eigs', verbosity=True)
v2 = np.real(g2.eigvect[0:slab.NT, 0])

# k3 = gamma(slab, PN, nev=nev, algo='eig', verbosity=True)
# v3 = np.real(k3.eigvect[0:slab.NT, 0])


fig = plt.figure()
plt.plot(slab.mesh, v1/norm(v1))
plt.plot(slab.mesh, v2/norm(v2), linestyle='dashed')
# plt.plot(slab.mesh, v3/norm(v3), linestyle='dotted')
plt.legend(('PETSc', 'scipy.sparse.eigs', 'scipy.linalg.eig'))

print(*1.8/g1.eigvals[0::2], sep='\n')
print(*1.8/g1.eigvals[0::2], sep='\n')
