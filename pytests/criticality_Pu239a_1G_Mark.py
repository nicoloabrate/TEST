"""
Created on Sat May 16 10:11:57 2020

author: N. Abrate.

file: .py

description:
"""
import os
import sys
sys.path.append('/home/caron/shared_memory/codici/python/TEST')
import time as t
from TEST.geometry import Slab
import TEST.NeutronTransportEquation as NTE
from TEST.eigenproblems.alpha import alphaprompt
from TEST.eigenproblems.kappa import kappa
from TEST.eigenproblems.gamma import gamma
from TEST.eigenproblems.delta import delta

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm

nev = 1
M = 50
G = 1
N = 11
bc = 'Mark'
H = 0.77032  # 0.76378  #


matname = ['Pu239a1L']
xlayers = [-H, H]
# define geometry and mesh
slab = Slab(M, xlayers, matname, [bc], G)
start = t.time()
PN = NTE.PN(slab, N, fmt='csc', allope=True)

print('Elapsed time PN: %f' % (t.time()-start))

k = kappa(slab, PN, nev=nev, algo='eigs', verbosity=True)
g = gamma(slab, PN, nev=nev+1, algo='eig', verbosity=True)
# d = delta(slab, PN, nev=nev, algo='eig', verbosity=True)
a = alphaprompt(slab, PN, nev=nev+5, algo='eigs', verbosity=True)

kv = np.real(k.eigvect[0:slab.NT, 0])
gv = np.real(g.eigvect[0:slab.NT, 0])
# dv = np.real(d.eigvect[0:slab.NT, 0])
av = np.real(a.eigvect[0:slab.NT, 0])


plt.rcParams.update({'font.size': 22})
fig = plt.figure()
plt.plot(slab.mesh, kv)
plt.plot(slab.mesh, gv, linestyle='dashed')
plt.plot(slab.mesh, av, linestyle='dotted')
# plt.plot(slab.mesh, dv, linestyle='dotted')
plt.legend(('$ \kappa -mode$', '$ \gamma -mode$', r'$ \alpha -mode$', '$\delta -mode$'))
plt.xlabel('x [cm]')
plt.ylabel('$\phi_0 $ [a.u.]')
plt.title(r'$ \alpha$=%.5e, $\kappa$=%.5f, $\gamma$=%.5f' % (a.eigvals[0], k.eigvals[0], g.eigvals[0]))
# plt.savefig('critical_modes_%s.pdf' % bc, bbox_inches = 'tight')
