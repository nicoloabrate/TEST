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

nev = 1
M = -10
N = 1
G = 1
bc = 'Mark'
H = 0.1
xlayers = [0, H]
# define geometry and mesh
myslab = Slab(M, xlayers, ['tmp2'], [bc], G, N, 'FD')
myslabD = Slab(M, xlayers, ['tmp2'], ['zero'], G, 0, 'FD')
myslab2 = Slab(M, xlayers, ['tmp2'], [bc], G, N+1, 'FD')

myDiff = NTE.Diffusion(myslabD, steady=True, fmt='csc', BC=True)
myPN = NTE.PN(myslab, N, steady=True, fmt='csc', BC=True)
mySN = NTE.SN(myslab2, N+1, steady=True, fmt='csc', BC=True)

# transform SN operators into PN operators to check consistency
S0 = 0.1
mysrc = lambda x: S0 # *(x>=H/4 and x<=3*H/4)
sD = sourceproblem(myDiff, 'static', myslabD, mysrc)
sP = sourceproblem(myPN, 'static', myslab, mysrc)
sS = sourceproblem(mySN, 'static', myslab2, mysrc)
sD.solve()
sP.solve()
sS.solve()

phiSN = sS.solution.flux
totphiSN = phiSN[0:myslab2.nS-1]+np.flipud(phiSN[myslab2.nS-1:2*myslab2.nS-1])

phi1 = lambda x: 20*(S0/(1*2)*(1-np.exp(-np.sqrt(3)*x)))
phi2 = lambda x: 20*(S0/(1*2)*(1-np.exp(-np.sqrt(3)*(H-x))))
# phi2 = lambda x: 20*(S0/(1*2)*(np.exp(np.sqrt(3)*x)-np.exp(np.sqrt(3)*(2*x-H))))
phiT = lambda x: phi1(x)+phi2(x)


plt.plot(myslab.mesh, np.flipud(phi1(myslab.mesh))+phi1(myslab.mesh), label='phitot')
# plt.plot(myslab.mesh, np.flipud(phiSN[myslab2.nS:2*myslab2.nS]))
plt.legend()
plt.show()

dx = myslab.dx[0]
phi1_num = np.zeros((myslab.nS, ))
phi2_num = np.zeros((myslab.nS, ))
q = (S0/2)*20
mu = 1/np.sqrt(3)
for i in range(1, myslab.nS):
    phi1_num[i] = q/(mu/dx+1/2)+phi1_num[i-1]*(1-dx/(2*mu))/(1+dx/(2*mu))
    phi2_num[myslab.nS-i-1] = q/(mu/dx+1/2)+phi2_num[myslab.nS-i]*(1-dx/(2*mu))/(1+dx/(2*mu))

# plt.plot(myslabD.mesh, sD.solution.flux, label='Diff')
plt.plot(myslab.mesh, sP.solution.flux[0:myslab.nS], label='PN-FD', c='r')
plt.plot(myslab2.stag_mesh, totphiSN, label='SN-FV')
plt.plot(myslab.mesh, phiT(myslab.mesh), '--', label='exact', c='b')
plt.plot(myslab.mesh, phi1_num+phi2_num, label='SN-FD')
# plt.plot(np.linspace(0, H/2, 1000), tmp, label='Diff. analytical')
plt.legend()
plt.show()
