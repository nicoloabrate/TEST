#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 13:00:55 2021

@author: nabrate
"""
import sys
sys.path.append('/opt/programs_nabrate/mycodes')
import numpy as np
from TEST.geometry import Slab
import TEST.models.NeutronTransportEquation as NTE
from TEST.models.EigenProblem import eigenproblem
from TEST.models.SourceProblem import sourceproblem
import matplotlib.pyplot as plt

H = 5
R = 0.5
xlayers = [-R, R] # -H, H,
M = -8 # [-5, -5, -5]
N = 1
G = 2
bc = 'Mark'

# define geometry and mesh
mats = ['ideal_reflector'] # ['H2O', 'H2O', 'H2O']
# myslabD = Slab(M, xlayers, mats, ['zero'], G, 0, 'FD')
myslabP = Slab(M, xlayers, mats, [bc], G, N, 'FD')
myslabS = Slab(M, xlayers, mats, [bc], G, N+1, 'FD')

# myDiff = NTE.Diffusion(myslabD, steady=True, fmt='csc', BC=True)
myPN = NTE.PN(myslabP, N, steady=True, fmt='csc', BC=True)
mySN = NTE.SN(myslabS, N+1, steady=True, fmt='csc', BC=True)

# transform SN operators into PN operators to check consistency
S0 = 1
mysrc = lambda x: S0 # , E: S0*(E>0.0625)+S0/2*(E<=0.0625) # *(x>=-H and x<=H)
# build problem loading the operators
# sD = sourceproblem(myDiff, 'static', myslabD, mysrc)
sP = sourceproblem(myPN, 'static', myslabP, mysrc)
sS = sourceproblem(mySN, 'static', myslabS, mysrc)
# solve problem
# sD.solve()
sP.solve()
sS.solve()

# plot solutions
# fig, ax = plt.subplots()
# myslabD.displaygeom()
# sD.solution.plot(1)
# if G > 1:
#     sD.solution.plot(2)
# plt.title('Diffusion')

fig, ax = plt.subplots()
myslabP.displaygeom()
sP.solution.plot(1)
if G > 1:
    sP.solution.plot(2)
plt.title('P1')

fig, ax = plt.subplots()
myslabS.displaygeom()
sS.solution.plot(1)
if G > 1:
    sS.solution.plot(2)
plt.title('S2')

# ---
fig, ax = plt.subplots()
# phiD = sD.solution.get(group=1, moment=0, mode=0)
phiP1 = sP.solution.get(group=1, moment=0, mode=0)
phiS2 = sS.solution.get(group=1, moment=0, mode=0)

# plt.plot(myslabD.mesh, phiD, label='Diffusion, g=1')
plt.plot(myslabP.mesh, phiP1, label='P1, g=1')
plt.plot(myslabS.mesh, phiS2, label='S2, g=1', linestyle='--')

if G > 1:
    # phiD = sD.solution.get(group=2, moment=0, mode=0)
    phiP1 = sP.solution.get(group=2, moment=0, mode=0)
    phiS2 = sS.solution.get(group=2, moment=0, mode=0)

    # plt.plot(myslabD.mesh, phiD, label='Diffusion, g=2')
    plt.plot(myslabP.mesh, phiP1, label='P1, g=2')
    plt.plot(myslabS.mesh, phiS2, label='S2, g=2', linestyle='--')

plt.legend()
plt.show()

# check maxima
# print('D  {:3f}'.format(phiD.max()))
print('P1 {:3f}'.format(phiP1.max()))
print('S2 {:3f}'.format(phiS2.max()))

# get operators
# Lp1=myPN.L.todense()
# Sp1=myPN.S.todense()
# Fp1=myPN.F.todense()
# Rp1=myPN.R.todense()
# S0p1=myPN.S0.todense()
# F0p1=myPN.F0.todense()

Ls2=mySN.L.todense()
Ss2=mySN.S.todense()
Fs2=mySN.F.todense()
Rs2=mySN.R.todense()
S0s2=mySN.S0.todense()
F0s2=mySN.F0.todense()

# plot angular flux
mu = myslabS.QW['mu']
fig, ax = plt.subplots()
# phiD = sD.solution.get(group=1, moment=0, mode=0)
phi_mu1 = sS.solution.get(group=1, angle=mu[0], mode=0)
phi_mu2 = sS.solution.get(group=1, angle=mu[1], mode=0)

# plt.plot(myslabD.mesh, phiD, label='Diffusion, g=1')
plt.plot(myslabS.mesh, phi_mu1, label='mu1, g=1')
plt.plot(myslabS.mesh, phi_mu2, label='mu2, g=1', linestyle='--')

if G > 1:
    phi_mu1 = sS.solution.get(group=2, angle=mu[0], mode=0)
    phi_mu2 = sS.solution.get(group=2, angle=mu[1], mode=0)

    plt.plot(myslabS.mesh, phi_mu1, label='mu1, g=2')
    plt.plot(myslabS.mesh, phi_mu2, label='mu2, g=2', linestyle='--')

plt.legend()
plt.show()