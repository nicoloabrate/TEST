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
myslabS = Slab(M, xlayers, ['tmp2'], [bc], G, N+1, 'FD')
myslabP = Slab(M, xlayers, ['tmp'], [bc], G, N, 'FD')

myPN = NTE.PN(myslabP, N, steady=True, fmt='csc', BC=True)
mySN = NTE.SN(myslabS, N+1, steady=True, fmt='csc', BC=True)

kP = eigenproblem(myPN, 'kappa', myslabP)
kP.solve(algo='eig')
Ak = kP.A

mysrc = np.ones(((myslabS.nS-1)*2, ))
dx = myslabS.dx[0]
for i in range(0, 3):
    
    # transform SN operators into PN operators to check consistency
    sS = sourceproblem(mySN, 'static', myslabS, mysrc)
    sS.solve()
    phiSN = sS.solution.flux
    
    phi1_num = np.zeros((myslabS.nS, ))
    phi2_num = np.zeros((myslabS.nS, ))
    q = mysrc
    mu = 1/np.sqrt(3)
    for i in range(1, myslabS.nS):
        phi1_num[i] = q[i]/(mu/dx+1/2)+phi1_num[i-1]*(1-dx/(2*mu))/(1+dx/(2*mu))
        phi2_num[myslabS.nS-i-1] = q[i]/(mu/dx+1/2)+phi2_num[myslabS.nS-i]*(1-dx/(2*mu))/(1+dx/(2*mu))

    
    # update source
    phiPN, _ = kP.solution.get(0)
    phi = np.interp(myslabS.stag_mesh, myslabS.mesh, phiPN)
    mysrc = np.concatenate([phi*1/2, phi*1/2])/kP.solution.eigvals.real[0]

    mysrcP = np.zeros((myslabS.nS*2-1, ))
    if i == 0:
        mysrcP[0:myslabS.nS] = 1
    else:
        mysrcP[0:myslabS.nS] = phiPN/kP.solution.eigvals.real[0]

    sP = sourceproblem(myPN, 'static', myslabP, mysrcP)
    sP.solve()
    As = sP.A

    totphiSN = phiSN[0:myslabS.nS-1]+np.flipud(phiSN[myslabS.nS-1:2*myslabS.nS-1])
    plt.plot(myslabP.mesh, kP.solution.eigvect[0:myslabS.nS, 0], label='PN-FD-k', c='r')
    plt.plot(myslabP.mesh, sP.solution.flux[0:myslabS.nS], label='PN-FD', c='g')
    plt.plot(myslabS.stag_mesh, totphiSN, label='SN-FV')
    plt.plot(myslabS.mesh, phi1_num+phi2_num, label='SN-FD')
    plt.legend()
    plt.show()

    # prova ad interpolare la sorgente dalla P1 per verificare che i profili si sovrappongano
    plt.plot(myslabS.stag_mesh, phiSN[0:myslabS.nS-1], label='phi1 SN-FV')
    plt.plot(myslabS.mesh, phi1_num[0:myslabS.nS], label='phi1 SN-FD')
    plt.legend()
    plt.show()