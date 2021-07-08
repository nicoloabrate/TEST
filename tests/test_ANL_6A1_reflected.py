"""
Author: N. Abrate.

File: test_kDiffusion_analyticalBenchmark_refl2G.py

Description: Analytical benchmark for a reflected 2G slab. Both forward and
            adjoint models are executed.
"""
import sys
sys.path.append('../../')
# import pytest
from TEST.geometry import Slab
import TEST.models.NeutronTransportEquation as NTE
import TEST.models.AdjointTransportEquation as ATE
from TEST.models.EigenProblem import eigenproblem
# import matplotlib.pyplot as plt
import time as t
import numpy as np

kref = 0.9015507
matfuel = 'ANL_6A1_R1_R3'
matrefl = 'ANL_6A1_R2'
G = 2
xlayers = [-120, -80, 80, 120]
mats = ['ANL_6A1_R1_R3', 'ANL_6A1_R2', 'ANL_6A1_R1_R3']

algo = 'eig'
nev = 1
M = [-21, -81, -21]
N = 0
bc = 'zero'

# P1
N = 0
bc = 'zero'
myslabD = Slab(M, xlayers, mats, [bc], G, N, 'FD')

t0 = t.time()
myDiff = NTE.Diffusion(myslabD, steady=True, fmt='csc')
kD = eigenproblem(myDiff, 'kappa', myslabD, nev=nev)
print('Elapsed time to setup PN FD: {} s'.format(t.time()-t0))

t0 = t.time()
kD.solve(algo=algo)
print('Elapsed time, diffusion: {} s'.format(t.time()-t0))
print('kD={:5f}'.format(kD.solution.eigvals[0]))

# P1
N = 1
bc = 'Mark'
myslabP = Slab(M, xlayers, mats, [bc], G, N, 'FD')

t0 = t.time()
myPN = NTE.PN(myslabP, N, steady=True, fmt='csc')
kP1 = eigenproblem(myPN, 'kappa', myslabP, nev=nev)
print('Elapsed time to setup PN FD: {} s'.format(t.time()-t0))

t0 = t.time()
kP1.solve(algo=algo)
print('Elapsed time, PN: {} s'.format(t.time()-t0))
print('kP1={:5f}'.format(kP1.solution.eigvals[0]))

# S2 FD
N = N+1
bc = 'Mark'
myslabSFD = Slab(M, xlayers, mats, [bc], G, N, 'FD')

t0 = t.time()
mySNFD = NTE.SN(myslabSFD, N, steady=True, fmt='csc', BC=True)
kS2FD = eigenproblem(mySNFD, 'kappa', myslabSFD, nev=nev)
print('Elapsed time to setup SN FD: {} s'.format(t.time()-t0))

t0 = t.time()
kS2FD.solve(algo=algo)
print('Elapsed time, SN FD: {} s'.format(t.time()-t0))
print('kS2 FD={:5f}'.format(kS2FD.solution.eigvals[0]))

# S2 FV
bc = 'Mark'
myslabSFV = Slab(M, xlayers, mats, [bc], G, N, 'FV')

t0 = t.time()
mySNFV = NTE.SN(myslabSFV, N, steady=True, fmt='csc', BC=True)
kS2FV = eigenproblem(mySNFV, 'kappa', myslabSFV, nev=nev)
print('Elapsed time to setup SN FV: {} s'.format(t.time()-t0))

t0 = t.time()
kS2FV.solve(algo=algo)
print('Elapsed time, SN FV: {} s'.format(t.time()-t0))
print('kS2 FV={:5f}'.format(kS2FV.solution.eigvals[0]))


phiP1_g1 = kP1.solution.get(group=1, moment=0, mode=0, normalise='peaktotalflux')
phiS2FD_g1 = kS2FD.solution.get(group=1, moment=0, mode=0, normalise='peaktotalflux')
phiS2FV_g1 = kS2FV.solution.get(group=1, moment=0, mode=0, normalise='peaktotalflux')

# get operators
Lp1=np.asarray(myPN.L.todense())
Sp1=np.asarray(myPN.S.todense())
Fp1=np.asarray(myPN.F.todense())
Rp1=np.asarray(myPN.R.todense())
S0p1=np.asarray(myPN.S0.todense())
F0p1=np.asarray(myPN.F0.todense())

L_FD=np.asarray(mySNFD.L.todense())
S_FD=np.asarray(mySNFD.S.todense())
F_FD=np.asarray(mySNFD.F.todense())
R_FD=np.asarray(mySNFD.R.todense())
S0_FD=np.asarray(mySNFD.S0.todense())
F0_FD=np.asarray(mySNFD.F0.todense())

L_FV=np.asarray(mySNFV.L.todense())
S_FV=np.asarray(mySNFV.S.todense())
F_FV=np.asarray(mySNFV.F.todense())
R_FV=np.asarray(mySNFV.R.todense())
S0_FV=np.asarray(mySNFV.S0.todense())
F0_FV=np.asarray(mySNFV.F0.todense())