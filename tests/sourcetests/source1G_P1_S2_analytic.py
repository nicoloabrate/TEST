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
from TEST.models.SourceProblem import sourceproblem
import matplotlib.pyplot as plt
from matplotlib import rcParams, colors
rcParams['figure.dpi'] = 200
plt.rc('grid', linestyle=":", color='black')

nev = 1
M = -50
N = 1
G = 1
bc = 'Mark'
H = 3
xlayers = [0, H]
# define geometry and mesh
myslab = Slab(M, xlayers, ['FullAbsorption'], [bc], G, N, 'FD')
myslabD = Slab(M, xlayers, ['FullAbsorption'], ['zero'], G, 0, 'FD')
myslab2 = Slab(M, xlayers, ['FullAbsorption'], [bc], G, N+1, 'FD')

myDiff = NTE.Diffusion(myslabD, steady=True, fmt='csc', BC=True)
myPN = NTE.PN(myslab, N, steady=True, fmt='csc', BC=True)
mySN = NTE.SN(myslab2, N+1, steady=True, fmt='csc', BC=True)

# transform SN operators into PN operators to check consistency
S0 = 100
mysrc = lambda x: S0
sD = sourceproblem(myDiff, 'static', myslabD, mysrc)
sP = sourceproblem(myPN, 'static', myslab, mysrc)
sS = sourceproblem(mySN, 'static', myslab2, mysrc)
sD.solve()
sP.solve()
sS.solve()

L = myslabD.regions['FullAbsorption'].getxs('DiffLength')
XSA = myslabD.regions['FullAbsorption'].getxs('Abs')
a = H/4
b = H/2
A = S0/XSA/(1/np.tanh((b-a)/L)*np.cosh(a/L)+np.sinh(a/L))
C = -A*np.cosh(a/L)/np.sinh((b-a)/L)
diff_anly = lambda x: (x>=0 and x<=H/4)*A*np.sinh(x/L)+(x>=H/4 and x<=H/2)*(S0/XSA+C*np.cosh((b-x)/L))
tmp = []
for x in np.linspace(0, H/2, 1000):
    tmp.append(diff_anly(x)[0][0]*20)

# source problem
phiSN = sS.solution.flux
totphiSN = phiSN[0:myslab2.nS]+np.flipud(phiSN[myslab2.nS:2*myslab2.nS])

phi1 = lambda x: 20*(S0/(XSA*2)*(1-np.exp(-np.sqrt(3)*XSA*x)))
phi2 = lambda x: 20*(S0/(XSA*2)*(1-np.exp(-np.sqrt(3)*XSA*(H-x))))
phiT = lambda x: phi1(x)+phi2(x)

plt.plot(myslab2.mesh, phiSN[0:myslab2.nS], label='phi1 SN')
plt.plot(myslab.mesh, phi1(myslab.mesh), label='phi1 exact')
plt.plot(myslab.mesh, np.flipud(phi1(myslab.mesh)), label='phi1 flipped')
plt.plot(myslab.mesh, phi2(myslab.mesh), label='phi2 exact')
plt.plot(myslab2.mesh, np.flipud(phiSN[myslab2.nS:]), label='phi2 SN')
plt.legend()
plt.show()

plt.plot(myslab.mesh, np.flipud(phi1(myslab.mesh))+phi1(myslab.mesh), label='phitot')
# plt.plot(myslab.mesh, np.flipud(phiSN[myslab2.nS:2*myslab2.nS]))
plt.legend()
plt.show()

# plt.plot(myslabD.mesh, sD.solution.flux, label='Diff')
plt.plot(myslab.mesh, sP.solution.flux[0:myslab.nS], label='$P_{1}$', c='r')
plt.plot(myslab2.mesh, totphiSN, label='$S_{2}$')
plt.plot(myslab.mesh, phiT(myslab.mesh), '--', label='analytic', c='b')
# plt.plot(np.linspace(0, H/2, 1000), tmp, label='Diff. analytical')
plt.ylabel('Flux [n/(cm$^2$ s)]')
plt.xlabel('x [cm]')
plt.legend()
plt.grid()
plt.savefig('source_analytic.png')
plt.show()

dx = myslab.dx[0]
phi1_num = np.zeros((myslab.nS, ))
phi2_num = np.zeros((myslab.nS, ))
q = (S0/2)*20
mu = 1/np.sqrt(3)
for i in range(1, myslab.nS):
    phi1_num[i] = q/(mu/dx+XSA/2)+phi1_num[i-1]*(mu/dx-XSA/2)/(mu/dx+XSA/2)
    phi2_num[myslab.nS-i-1] = q/(mu/dx+XSA/2)+phi2_num[myslab.nS-i]*(mu/dx-XSA/2)/(mu/dx+XSA/2)


plt.plot(myslab2.mesh, phiSN[0:myslab2.nS], label='phi1 SN')
plt.plot(myslab.mesh, phi1(myslab.mesh), label='phi1 exact')
plt.plot(myslab.mesh, phi1_num, label='phi1 SN+FD')
plt.legend()
plt.show()

plt.plot(myslab2.mesh, np.flipud(phiSN[myslab2.nS:2*myslab2.nS]), label='phi2 SN')
plt.plot(myslab.mesh, phi2(myslab.mesh), label='phi2 exact')
plt.plot(myslab.mesh, phi2_num, label='phi2 SN+FD')
plt.legend()
plt.show()