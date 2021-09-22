#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 13:00:55 2021

@author: nabrate
"""
import sys
sys.path.append('/opt/programs_nabrate/mycodes')
sys.path.append('C:\\Users\\39346\\Documenti\\mycodes\\TEST')
from TEST.geometry import Slab
import TEST.models.NeutronTransportEquation as NTE
from TEST.models.EigenProblem import eigenproblem
import matplotlib.pyplot as plt

algo = 'eig'
M = [-4, -7]
N = 1
G = 1
bc = 'Mark'

# # this works
# mat2 = 'Absorber_0L'
# mat1 = 'Pu239_0L_noS'
# H = 2
# b = 1
# xlayers = [-b, b, H]

# this does not work
mat1 = 'Absorber_0L' # 'Pu239_0L_noFS' #
mat2 = 'Pu239_0L_noS' # 'MontagniniFuel' #
H = 2
b = 1
xlayers = [-b, b, H]

# mats = [mat2, mat1, mat2]
mats = [mat1, mat2]
# define geometry and mesh
myslabP = Slab(M, xlayers, mats, [bc], G, N, 'FD')
myslabSFV = Slab(M, xlayers, mats, [bc], G, N+1, 'FV')  #[mat1, mat2, mat1]
myslabSFD = Slab(M, xlayers, mats, [bc], G, N+1, 'FD')

myPN = NTE.PN(myslabP, N, steady=True, fmt='csc', BC=True)
mySNFV = NTE.SN(myslabSFV, N+1, steady=True, fmt='csc', BC=True)
mySNFD = NTE.SN(myslabSFD, N+1, steady=True, fmt='csc', BC=True)

# transform SN operators into PN operators to check consistency
sP = eigenproblem(myPN, 'kappa', myslabP)
sS_FV = eigenproblem(mySNFV, 'kappa', myslabSFV)
sS_FD = eigenproblem(mySNFD, 'kappa', myslabSFD)

sP.solve(normalisation='peaktotalflux', algo=algo)
sS_FV.solve(normalisation='peaktotalflux', algo=algo)
sS_FD.solve(normalisation='peaktotalflux', algo=algo)

print('Fundamental k-eigenvalue:')
print('P{}: {}'.format(N, sP.solution.eigvals[0]))
print('S{} FV: {}'.format(N+1, sS_FV.solution.eigvals[0]))
print('diff: {} [pcm]'.format(1E5*abs(sS_FV.solution.eigvals[0]-sP.solution.eigvals[0]).real))

print('S{} FD: {}'.format(N+1, sS_FD.solution.eigvals[0]))
print('diff: {} [pcm]'.format(1E5*abs(sS_FD.solution.eigvals[0]-sP.solution.eigvals[0]).real))

phiSN1_FV = sS_FV.solution.get(group=1, moment=0, mode=0)
phiSN1_FD = sS_FD.solution.get(group=1, moment=0, mode=0)
phiPN1 = sP.solution.get(group=1, moment=0, mode=0)

# plt.plot(myslabD.mesh, sD.solution.eigvect[0:myslabS.nS, 0], label='Diffusion')
plt.plot(myslabP.mesh, phiPN1, label='P1')
plt.plot(myslabSFV.mesh, phiSN1_FV, label='SN-FV', linestyle='--')
plt.plot(myslabSFD.mesh, phiSN1_FD, label='SN-FD', linestyle=':')
plt.legend()
plt.show()


L_FD = mySNFD.L.todense()
S_FD = mySNFD.S.todense()
F_FD = mySNFD.F.todense()
R_FD = mySNFD.R.todense()
S0_FD = mySNFD.S0.todense()
F0_FD = mySNFD.F0.todense()

L_FV = mySNFV.L.todense()
S_FV = mySNFV.S.todense()
F_FV = mySNFV.F.todense()
R_FV = mySNFV.R.todense()
S0_FV = mySNFV.S0.todense()
F0_FV = mySNFV.F0.todense()