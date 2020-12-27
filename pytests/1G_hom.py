"""
Created on Sat May 16 10:11:57 2020

author: N. Abrate.

file: .py

description:
"""
from test.geometry import Slab
import test.NeutronTransportEquation as NTE
import test.eigenproblems.kappa as kappa
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm

nev = 1
M = 10
G = 1
N = 3
bc = 'Mark'
H = 10

matname = ['Pu239']
xlayers = [0, H]
# define geometry and mesh
slab = Slab(M, xlayers, matname, [bc], G)
PN = NTE.PN(slab, N, fmt='csc')

k1 = kappa(slab, PN, nev=nev, verbosity=True)
v1 = np.real(k1.eigvect[0:slab.NT, 0])

k2 = kappa(slab, PN, nev=nev, algo='eigs', verbosity=True)
v2 = np.real(k2.eigvect[0:slab.NT, 0])

# k3 = kappa(slab, PN, nev=nev, algo='eig', verbosity=True)
# v3 = np.real(k3.eigvect[0:slab.NT, 0])


fig = plt.figure()
plt.plot(slab.mesh, v1/norm(v1))
plt.plot(slab.mesh, v2/norm(v2), linestyle='dashed')
# plt.plot(slab.mesh, v3/norm(v3), linestyle='dotted')
plt.legend(('PETSc', 'scipy.sparse.eigs', 'scipy.linalg.eig'))