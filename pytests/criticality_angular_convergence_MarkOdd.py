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
from test.geometry import Slab
import test.NeutronTransportEquation as NTE
from test.eigenproblems.alpha import alphaprompt
from test.eigenproblems.kappa import kappa
from test.eigenproblems.gamma import gamma
from test.utils import h5
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm

nev = 1
M = 100
G = 1
bc = 'MarkOdd'
H = 0.77032
Nmax = 40

k_ref = 1.0000039336446056
g_ref = 1.0000019020968423
a_ref = -1521.6736698184789

matname = ['Pu239a1L']
xlayers = [-H, H]
# define geometry and mesh
slab = Slab(M, xlayers, matname, [bc], G)
start = t.time()

kappa0 = []
gamma0 = []
alpha0 = []
for n in range(1, Nmax):
    N = n
    PN = NTE.PN(slab, N, fmt='csc', allope=True)
    k = kappa(slab, PN, nev=nev, verbosity=True)
    g = gamma(slab, PN, nev=nev, verbosity=True)
    a = alphaprompt(slab, PN, nev=nev, algo='eig', verbosity=True)
    kappa0.append(k.eigvals[0])
    gamma0.append(g.eigvals[0])
    a0 = np.min(abs(a.eigvals[a.eigvals.imag == 0].real))
    pos = np.where(abs(a.eigvals) == a0)
    alpha0.append(np.asscalar(a.eigvals[abs(a.eigvals) == a0, ]))


dk = 1e5*(np.asarray(kappa0)-k_ref).real
dg = 1e5*(np.asarray(gamma0)-g_ref).real
da = (np.asarray(alpha0)-a_ref).real

handles = []
fig = plt.figure()
ax = plt.gca()
ax.yaxis.grid(True)
b1 = plt.bar(range(1, Nmax, 2), dk[0::2], color='midnightblue',
             width=1, edgecolor='black', align='center')
handles.append(b1)
b2 = plt.bar(range(2, Nmax-1, 2), dk[1::2], color='seagreen',
             width=1, edgecolor='black', align='center')
handles.append(b2)
plt.yscale('symlog')
plt.legend((b1, b2), ('odd $P_N$', 'even $P_N$'))
plt.xlim((0.5, Nmax-0.5))
plt.ylabel('$\kappa_{N}-\kappa_{1001}$ [pcm]')
labels_odd = [str(val) for val in range(1, Nmax, 2)]
plt.xticks(ticks=range(1, Nmax, 2), labels=labels_odd)
ax.set_yticks([-1e5, -1e3, -1e1, 0, 1e1, 1e3, 1e5])
plt.savefig('hist_kappa_%s_%g.pdf' % (bc, M))

handles = []
fig = plt.figure()
ax = plt.gca()
ax.yaxis.grid(True)
b1 = plt.bar(range(1, Nmax, 2), dg[0::2], color='midnightblue',
             width=1, edgecolor='black', align='center')
handles.append(b1)
b2 = plt.bar(range(2, Nmax-1, 2), dg[1::2], color='seagreen',
             width=1, edgecolor='black', align='center')
handles.append(b2)
plt.yscale('symlog')
plt.legend((b1, b2), ('odd $P_N$', 'even $P_N$'))
plt.xlim((0.5, Nmax-0.5))
plt.ylabel('$\gamma_{N}-\gamma_{1001}$ [pcm]')
labels_odd = [str(val) for val in range(1, Nmax, 2)]
plt.xticks(ticks=range(1, Nmax, 2), labels=labels_odd)
ax.set_yticks([-1e5, -1e3, -1e1, 0, 1e1, 1e3, 1e5])
plt.savefig('hist_gamma_%s_%g.pdf' % (bc, M))

handles = []
fig = plt.figure()
ax = plt.gca()
ax.yaxis.grid(True)
b1 = plt.bar(range(1, Nmax, 2), da[0::2].T, color='midnightblue',
             width=1, edgecolor='black', align='center')
handles.append(b1)
b2 = plt.bar(range(2, Nmax-1, 2), da[1::2], color='seagreen',
             width=1, edgecolor='black', align='center')
handles.append(b2)
plt.yscale('symlog')
plt.legend((b1, b2), ('odd $P_N$', 'even $P_N$'))
plt.xlim((0.5, Nmax-0.5))
plt.ylabel(r'$ \alpha_{N}-\alpha_{1001} [s^{-1}]$')
labels_odd = [str(val) for val in range(1, Nmax, 2)]
plt.xticks(ticks=range(1, Nmax, 2), labels=labels_odd)
ax.set_yticks([-1e5, -1e3, -1e1, 0, 1e1, 1e3, 1e5])
plt.savefig('hist_alpha_%s_%g.pdf' % (bc, M))
