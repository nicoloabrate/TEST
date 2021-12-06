# %%
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
from scipy.special import roots_legendre, eval_legendre
from matplotlib import rc, rcParams
rcParams['figure.dpi'] = 100

M = -10
N = 2
G = 1
bc = 'Mark'
H = 0.1
xlayers = [0, H]
material = 'Absorber_0L'
scheme = 'FD'

# ensure positive and then negative directions
mu, w = roots_legendre(3)
mu[::-1].sort()

# define geometry and mesh
myslab = Slab(M, xlayers, [material], bc, G, N, 'FD')
myslabD = Slab(M, xlayers, [material], 'zero', G, 0, 'FD')
myslab2 = Slab(M, xlayers, [material], bc, G, N+1, scheme)

myDiff = NTE.Diffusion(myslabD, steady=True, fmt='csc', BC=True)
myPN = NTE.PN(myslab, N, steady=True, fmt='csc', BC=True)
mySN = NTE.SN(myslab2, N+1, steady=True, fmt='csc', BC=True)

# transform SN operators into PN operators to check consistency
S0 = 0.1
mysrc = lambda x: S0
q = S0*20
# sD = sourceproblem(myDiff, 'static', myslabD, mysrc)
sS = sourceproblem(mySN, 'static', myslab2, mysrc)
sP = sourceproblem(myPN, 'static', myslab, mysrc)
# sD.solve()
sS.solve()
sP.solve()

phiSN = sS.solution
# --- analytical
XS = myslab2.regions[material].Tot
MFP = 1/XS
phi1 = lambda x: q/2/XS*(1-np.exp(-XS*x/mu[0]))
phi2 = lambda x: q/2/XS*np.ones((len(x)))
phi3 = lambda x: q/2/XS*(1-np.exp(-XS*(H-x)/mu[0]))
xx = myslab.mesh
phi_tot3 = w[0]*phi1(xx)+w[1]*phi2(xx)+w[2]*phi3(xx)
# --- numerical
dx = myslab.dx[0]
phi1_num = np.zeros((myslab.nS, ))
phi2_num = np.zeros((myslab.nS, ))
phi3_num = np.zeros((myslab.nS, ))
phi2_num[0] = q/2/XS
for i in range(1, myslab.nS):
    phi1_num[i] = q/2/(mu[0]/dx+XS/2)+phi1_num[i-1]*(1-XS*dx/(2*mu[0]))/(1+XS*dx/(2*mu[0]))
    # take absolute value of mu (so mu0)
    phi3_num[myslab.nS-i-1] = q/2/(mu[0]/dx+XS/2)+phi3_num[myslab.nS-i]*(1-XS*dx/(2*mu[0]))/(1+XS*dx/(2*mu[0]))
    phi2_num[i] = q/2/XS

phi_tot_num = w[0]*phi1_num+w[1]*phi2_num+w[2]*phi3_num
# plot angular flux
plt.plot(xx, phi1_num, label='numerical $\mu_1$')
plt.plot(xx, phi2_num, label='numerical $\mu_2$')
plt.plot(xx, phi3_num, label='numerical $\mu_3$')
plt.plot(xx, phi1(xx)[0], ls='--', label='analytical $\mu_1$')
plt.plot(xx, phi2(xx)[0], ls='--', label='analytical $\mu_2$')
plt.plot(xx, phi3(xx)[0], ls='--', label='analytical $\mu_3$')
phiSN.xplot(1, angle=mu[0], ls='-.', label='TEST $\mu_1$')
phiSN.xplot(1, angle=mu[1], ls='-.', label='TEST $\mu_2$')
phiSN.xplot(1, angle=mu[2], ls='-.', label='TEST $\mu_3$')
plt.legend(bbox_to_anchor=(1.01, 1))
plt.show()


# plt.plot(myslabD.mesh, sD.solution.flux, label='Diff')
phiSN.xplot(1, label='TEST', lw=2)
plt.plot(xx, phi_tot_num, label='numerical', ls='--', lw=1)
plt.plot(xx, phi_tot3[0], label='analytical', ls='-.', lw=1)
plt.legend(bbox_to_anchor=(1.01, 1))
plt.show()

# # spy operators
# mySN.spy('C'); plt.title('C'); plt.show()
# mySN.spy('S0'); plt.title('S0'); plt.show()
# mySN.spy('F0'); plt.title('F0'); plt.show()
# mySN.spy('S'); plt.title('S'); plt.show()
# mySN.spy('F'); plt.title('F'); plt.show()
# mySN.spy('R'); plt.title('R'); plt.show()
# mySN.spy('L'); plt.title('L'); plt.show()

# %%
