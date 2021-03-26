#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 09:13:24 2021

@author: nabrate
"""
import sys
sys.path.append('../../')
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

M = 20
G = 2
N = 0
bc = 'zero'
H, R = 20, 40

nev = 10
matname = ['MontagniniFuel']
xlayers = [-R, R]
# define geometry and mesh
slab0 = Slab(M, xlayers, matname, [bc], G, 0, 'FD')
Fwd = NTE.Diffusion(slab0, steady=True, fmt='csc')
Adj = ATE.Diffusion(slab0, steady=True, fmt='csc')
# solve forward and adjoint problems
forward = kappa(slab0, Fwd, nev=nev, verbosity=True)
adjoint = kappa(slab0, Adj, nev=nev, verbosity=True)

# --- perturbation
slab = Slab(M, xlayers, matname, [bc], G, 0, 'FD')
perturbation = {'Fiss': {'where': [(-10, 0)], 'howmuch': [0, 3E-2]}}
slab.perturb(perturbation)
slab.displaygeom()
# --- solve perturbed problem
DiffP = NTE.Diffusion(slab, steady=True, fmt='csc')
pert = kappa(slab, DiffP, nev=nev, verbosity=True)

# --- plot
fig, ax = plt.subplots()
forward.plot(slab0, 1)
forward.plot(slab0, 2)
pert.plot(slab, 1, ls='--')
pert.plot(slab, 2, ls='--')
plt.legend(('$\phi_{ref, g=1}$', '$\phi_{ref, g=2}$', r'$\phi_{pert, g=1}$',
            r'$\phi_{pert, g=1}$'), loc='upper left', bbox_to_anchor=(1.05, 1))
plt.xlabel('x [cm]')
plt.ylabel('[a.u.]')
plt.ylim([0, 0.045])

# --- GPT
# GPT(4, slab, unperturbed, slab0, perturbed, forward, adjoint)