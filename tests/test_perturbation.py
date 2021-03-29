#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 09:13:24 2021

@author: nabrate
"""
import sys
sys.path.append('C:\\Users\\39346\\Documenti\\mycodes')
import numpy as np
import pytest
from TEST.geometry import Slab
import TEST.NeutronTransportEquation as NTE
import TEST.AdjointTransportEquation as ATE
from TEST.eigenproblems.criticality import kappa
from TEST.methods.GPT import GPT

import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
rcParams['figure.dpi'] = 200
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

delta = 0.1
M = 40
G = 2
N = 0
bc = 'zero'
H, R = 20, 40

nev = 20
matname = ['MontagniniFuel']
xlayers = [-R, R]
# define geometry and mesh
slab0 = Slab(M, xlayers, matname, [bc], G, 0, 'FD')
Fwd = NTE.Diffusion(slab0, steady=True, fmt='csc')
Adj = ATE.Diffusion(slab0, steady=True, fmt='csc')
# solve forward and adjoint problems
forward = kappa(slab0, Fwd, nev=nev)
forward.solve()
adjoint = kappa(slab0, Adj, nev=nev)
adjoint.solve()

# --- perturbation
slab = Slab(M, xlayers, matname, [bc], G, 0, 'FD')
perturbation = {'Nubar': {'where': [(-R, R)], 'howmuch': [0, delta]}}
slab.perturb(perturbation)
slab.displaygeom()
# --- solve perturbed problem
DiffP = NTE.Diffusion(slab, steady=True, fmt='csc')
pert = kappa(slab, DiffP, nev=nev)
pert.solve()

# --- plot
# fig, ax = plt.subplots()
# forward.plot(slab0, 1)
# forward.plot(slab0, 2)
# pert.plot(slab, 1, ls='--')
# pert.plot(slab, 2, ls='--')
# plt.legend(('$\phi_{ref, g=1}$', '$\phi_{ref, g=2}$', r'$\phi_{pert, g=1}$',
#             r'$\phi_{pert, g=1}$'), loc='upper left', bbox_to_anchor=(1.05, 1))
# plt.xlabel('x [cm]')
# plt.ylabel('[a.u.]')
# plt.ylim([0, 0.045])
m = slab.regions['Perturbation1']
L1 = m.DiffLength[0]
L2 = m.DiffLength[1]
B = np.pi/(2*R)
keff = m.Nsf[1]*m.S0[0, 1]/(m.Remxs[0]*m.Remxs[1]*(1+L1**2*B**2)*(1+L2**2*B**2))
# --- GPT
N = 10
gpt = GPT(N, slab0, forward, slab, pert, adjoint)

deltasum = [(-1)**n*delta**n for n in range(0, N)]
kp_ref = 1/(1/forward.eigvals[0]*sum(deltasum))
print((pert.eigvals[0]-gpt.pertEigv)*1E5)
print((kp_ref-gpt.pertEigv)*1E5)
print((kp_ref-pert.eigvals[0])*1E5)