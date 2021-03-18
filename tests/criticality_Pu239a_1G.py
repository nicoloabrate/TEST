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
from TEST.eigenproblems.kappa import kappa
from TEST.eigenproblems.gamma import gamma

import matplotlib.pyplot as plt
import numpy as np

nev = 1
M = 100
G = 1
N = 1
bc = 'Marshak'
H = 0.77032  # 0.76378


matname = ['Modak']
xlayers = [-H, H]
# define geometry and mesh
slab = Slab(M, xlayers, matname, [bc], G)
start = t.time()
PN = NTE.PN(slab, N, fmt='csc', allope=True)

print('Elapsed time PN: %f' % (t.time()-start))

k2 = kappa(slab, PN, nev=nev, verbosity=True)
g2 = gamma(slab, PN, nev=nev, verbosity=True)
# a1 = alphaprompt(slab, PN, nev=nev, verbosity=True, which='LR')
kv2 = np.real(k2.eigvect[0:slab.NT, 0])
gv2 = np.real(g2.eigvect[0:slab.NT, 0])
# av1 = np.real(a1.eigvect[0:slab.NT, 0])


k_eigs = kappa(slab, PN, nev=nev, algo='eigs', verbosity=True)
g_eigs = gamma(slab, PN, nev=nev, algo='eigs', verbosity=True)
# a2 = alphaprompt(slab, PN, nev=nev, algo='eigs', verbosity=True)
# kv2 = np.real(k2.eigvect[0:slab.NT, 0])
# gv2 = np.real(g2.eigvect[0:slab.NT, 0])
# av2 = np.real(a2.eigvect[0:slab.NT, 0])


# fig = plt.figure()
# plt.plot(slab.mesh, kv2)
# plt.plot(slab.mesh, gv2, linestyle='dashed')
# plt.plot(slab.mesh, av2, linestyle='dotted')
# plt.legend(('$ \kappa -mode$', '$ \gamma -mode$', r'$ \alpha -mode$'))
# plt.xlabel('x [cm]')
# plt.ylabel('$\phi_0 $ [a.u.]')
# plt.title(r'$ \alpha$=%.5e, $\kappa$=%.5f, $\gamma$=%.5f' % (a2.eigvals[0], k2.eigvals[0], g2.eigvals[0]))
# plt.savefig('subcritical_modes_%s.pdf' % bc)
