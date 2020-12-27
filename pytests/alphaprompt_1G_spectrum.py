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

nev = 1
M = 100
G = 1
bc = 'Mark'
H = 2
N = 7

matname = ['Pu239a1LforwhighS']
xlayers = [-H, H]
# define geometry and mesh
slab = Slab(M, xlayers, matname, [bc], G)
start = t.time()
PN = NTE.PN(slab, N, fmt='csc', steady=False, prompt=True)
print('Elapsed time PN: %f' % (t.time()-start))

a1 = alphaprompt(slab, PN, nev=nev, algo='eig', verbosity=True)
alpha = a1.eigvals[:]

a0 = np.min(abs(alpha[alpha.imag == 0].real))
pos = np.where(abs(alpha) == a0)
alpha0 = np.asscalar(alpha[abs(alpha) == a0, ])
alpha = np.delete(alpha, [pos])

fig = plt.figure()
plt.axhline(0, color='black', ls='dashed', lw = 0.5)
ax = plt.gca()
ax.scatter(alpha.real, alpha.imag, color='red')
ax.scatter(alpha0.real, alpha0.imag, color='blue', marker='*')
plt.ylabel(r'$Im(\alpha/v) ~[cm^{-1}]$')
plt.xlabel(r'$Re(\alpha/v) ~[cm^{-1}]$')
plt.ylim(-300, 300)
plt.xlim(-9.5, 0)
plt.savefig('spectrum_alpha_%s_%s_%g_%g_linlin.pdf' % (matname[0], bc, N, M))

fig = plt.figure()
plt.axhline(0, color='black', ls='dashed', lw = 0.5)
ax = plt.gca()
ax.scatter(alpha.real, alpha.imag, color='red')
ax.scatter(alpha0.real, alpha0.imag, color='blue', marker='*')
# plt.xscale('symlog')
plt.yscale('symlog')
plt.ylabel(r'$Im(\alpha/v) ~[cm^{-1}]$')
plt.xlabel(r'$Re(\alpha/v) ~[cm^{-1}]$')
ax.set_yticks([-1e3, -1e1, 0, 1e1, 1e3])
plt.ylim(-1e3, 1e3)
plt.xlim(-9.5, 0)
plt.savefig('spectrum_alpha_%s_%s_%g_%g_loglin.pdf' % (matname[0], bc, N, M))
