"""
Created on Sat May 16 10:11:57 2020

author: N. Abrate.

file: .py

description:
"""
import time as t
from test.geometry import Slab
import test.NeutronTransportEquation as NTE
from test.eigenproblems.alpha import alphaprompt
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm

nev = 1
M = 100
G = 1
bc = 'Mark'
H = 1
N = 15

k_ref = 1.0000039336446056
g_ref = 1.0000019020968423
a_ref = -1521.6736698184789

matname = ['Pu239a1L']
xlayers = [-H, H]
# define geometry and mesh
slab = Slab(M, xlayers, matname, [bc], G)
start = t.time()
PN = NTE.PN(slab, N, fmt='csc', steady=False, prompt=True)
print('Elapsed time PN: %f' % (t.time()-start))

a1 = alphaprompt(slab, PN, nev=nev, algo='eig', verbosity=True)
alpha = a1.eigvals[:]

fig = plt.figure()
ax = plt.gca()
ax.scatter(alpha.real, alpha.imag, color='red')
plt.yscale('symlog')
#plt.xscale('symlog')
plt.ylabel(r'$Im(\alpha_p) ~[s^{-1}]$')
plt.xlabel(r'$Re(\alpha_p) ~[s^{-1}]$')
ax.set_yticks([-1e7, -1e5, -1e3, -1e1, 0, 1e1, 1e3, 1e5, 1e7])
plt.ylim((-1e3, 1e3))
plt.savefig('spectrum_alpha_%s_%s_%g.pdf' % (matname[0], bc, M))