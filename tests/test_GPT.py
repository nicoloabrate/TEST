#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 09:13:24 2021

@author: nabrate
"""
import sys
sys.path.append('../')
sys.path.append('../../')
import numpy as np
import pytest
from TEST.geometry import Slab
import TEST.models.NeutronTransportEquation as NTE
import TEST.models.AdjointTransportEquation as ATE
from TEST.models.EigenProblem import eigenproblem
from TEST.methods.GPT import GPT
import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
rcParams['figure.dpi'] = 200
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
rcParams.update({'font.size': 18})

@pytest.mark.parametrize("delta", [(1E-3), (5E-3), (1E-2)])
def test_GPT_vs_analytical(delta):

    PO = 15
    nev = 10
    M = 100
    G = 2
    N = 0
    bc = 'zero'
    H = 40

    matname = ['MontagniniFuel']
    xlayers = [-H, H]
    # define geometry and mesh
    slab0 = Slab(M, xlayers, matname, [bc], G, N, 'FD')
    Fwd = NTE.Diffusion(slab0, steady=True, fmt='csc')
    Adj = ATE.Diffusion(slab0, steady=True, fmt='csc')
    # solve forward and adjoint problems
    forward = eigenproblem(Fwd, 'kappa', slab0, nev=nev)
    forward.solve()
    adjoint = eigenproblem(Adj, 'kappa', slab0, nev=nev)
    adjoint.solve()

    # --- perturbation
    slab = Slab(M, xlayers, matname, [bc], G, 0, 'FD')
    perturbation = {'Nubar': {'where': [(-H, H)], 'howmuch': [0, delta]}}
    slab.perturb(perturbation)
    slab.displaygeom()
    # --- initialise perturbed problem
    DiffP = NTE.Diffusion(slab, steady=True, fmt='csc')
    pert = eigenproblem(DiffP, 'kappa', slab, nev=nev)

    m = slab.regions['Perturbation1']
    L1 = m.DiffLength[0]
    L2 = m.DiffLength[1]
    B = np.pi/(2*H)
    # analytical perturbed eigenvalue
    keff = m.Nsf[1]*m.S0[0, 1]/(m.Remxs[0]*m.Remxs[1]*(1+L1**2*B**2)*(1+L2**2*B**2))
    # --- GPT
    gpt = GPT(PO, forward, pert, adjoint)
    # --- analytical perturbation
    deltasum = [(-1)**n*delta**n for n in range(0, PO)]
    kp_ref = 1/(1/forward.solution.eigvals[0]*sum(deltasum))

    assert abs((kp_ref-gpt.solution.eigvals[0])*1E5) < 10
    assert abs((keff-gpt.solution.eigvals[0])*1E5) < 10

@pytest.mark.parametrize("delta", [(1E-3), (5E-3), (1E-2)])
def test_GPT_vs_DirectPerturbation(delta):
    PO = 15
    nev = 15
    M = 20
    G = 2
    bc = 'zero'
    H = 40

    matname = ['MontagniniFuel']
    xlayers = [-H, H]
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
    perturbation = {'Nubar': {'where': [(0, H/2)], 'howmuch': [0, delta]}}
    slab.perturb(perturbation)
    slab.displaygeom()
    # --- solve perturbed problem
    DiffP = NTE.Diffusion(slab, steady=True, fmt='csc')
    pert = eigenproblem(DiffP, 'kappa', slab, nev=nev)
    pert.solve()

    # --- GPT
    gpt = GPT(PO, forward, pert, adjoint)

    # --- plot
    fig, ax = plt.subplots()
    gpt.solution.plot(1)
    gpt.solution.plot(2)
    pert.solution.plot(1, ls='--')
    pert.solution.plot(2, ls='--')
    plt.legend(('r$\phi_{GPT, g=1}$', 'r$\phi_{GPT, g=2}$', r'$\phi_{pert, g=1}$',
                r'$\phi_{pert, g=1}$'), loc='upper left', bbox_to_anchor=(1.05, 1))
    plt.xlabel('x [cm]')
    plt.ylabel('[a.u.]')

    assert (pert.solution.eigvals[0]-gpt.solution.eigvals[0])*1E5 < 10
