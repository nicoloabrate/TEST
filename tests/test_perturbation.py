#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 09:13:24 2021

@author: nabrate
"""
import sys
sys.path.append('../../')
import numpy as np
import pytest
from TEST.geometry import Slab
import TEST.NeutronTransportEquation as NTE
import TEST.AdjointTransportEquation as ATE
from TEST.eigenproblems.EigenProblem import eigenproblem
from TEST.methods.GPT import GPT

import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
rcParams['figure.dpi'] = 200
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

H = 8
nev = 1
M = 100
G = 1
N = 3
bc = 'Mark'
xlayers = [0, H]
# define geometry and mesh
myslab = Slab(M, xlayers, ['Modak'], [bc], G, N, 'FD')
myPN = NTE.PN(myslab, N, steady=True, fmt='csc')
k1 = eigenproblem(myPN, 'kappa', myslab, nev=nev)
k1.solve(verbosity=True, algo='eigs')

PO = 10
nev = 10
delta = 5E-2
M = 10
G = 2
N = 0
bc = 'zero'
H, R = 20, 40

matname = ['MontagniniFuel']
xlayers = [-R, R]
# define geometry and mesh
slab0 = Slab(M, xlayers, matname, [bc], G, 0, 'FD')
Fwd = NTE.Diffusion(slab0, steady=True, fmt='csc')
Adj = ATE.Diffusion(slab0, steady=True, fmt='csc')
# solve forward and adjoint problems
forward = eigenproblem(Fwd, 'kappa', slab0, nev=nev)
forward.solve()
adjoint = eigenproblem(Adj, 'kappa', slab0, nev=nev)
adjoint.solve()

# --- perturbation
slab = Slab(M, xlayers, matname, [bc], G, 0, 'FD')
perturbation = {'Nubar': {'where': [(0, R/2)], 'howmuch': [0, delta]}}
slab.perturb(perturbation)
slab.displaygeom()
# --- solve perturbed problem
DiffP = NTE.Diffusion(slab, steady=True, fmt='csc')
pert = eigenproblem(DiffP, 'kappa', slab, nev=nev)
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
gpt = GPT(PO, forward, pert, adjoint)

deltasum = [(-1)**n*delta**n for n in range(0, PO)]
kp_ref = 1/(1/forward.solution.eigvals[0]*sum(deltasum))
print((pert.solution.eigvals[0]-gpt.solution.eigvals[0])*1E5)
# print((kp_ref-gpt.pertEigv)*1E5)
# print((kp_ref-pert.eigvals[0])*1E5)