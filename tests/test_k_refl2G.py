"""
Author: N. Abrate.

File: test_kDiffusion_analyticalBenchmark_refl2G.py

Description: Analytical benchmark for a reflected 2G slab. Both forward and
            adjoint models are executed.
"""
import sys
sys.path.append('../../')
import pytest
from TEST.geometry import Slab
import TEST.models.NeutronTransportEquation as NTE
import TEST.models.AdjointTransportEquation as ATE
from TEST.models.EigenProblem import eigenproblem
import matplotlib.pyplot as plt


H, R, G, matrefl, ref = (5, 10, 1, 'Absorber_0L', 1.004241348107076)  # Absorber_0L
matfuel = 'Pu239a' #  'MontagniniFuel' # 'Pu239_0L_noS'# 'MontagniniFuel' # 'Pu239a' # 
G = 1
xlayers = [-R, R]
mats = [matfuel] # [matrefl, matfuel, matrefl, matfuel, matrefl] # [matfuel]*(len(xlayers)-1) # 
#(40, 100, 2, 'MontagniniReflector2', 1.045766651960480),
# (40, 100, 2, 'MontagniniReflector3', 1.020903109926185)])
algo = 'eig'
nev = 1
M = [-50]*(len(xlayers)-1)
N = 0
bc = 'zero'
# Diffusion
myslabD = Slab(M, xlayers, mats, [bc], G, N, 'FD')
myPN = NTE.Diffusion(myslabD, steady=True, fmt='csc')
kD = eigenproblem(myPN, 'kappa', myslabD, nev=nev)
kD.solve(algo=algo)
print('kD={:5f}'.format(kD.solution.eigvals[0]))

# P1
N = 1
bc = 'Mark'
myslabP = Slab(M, xlayers, mats, [bc], G, N, 'FD')
myPN = NTE.PN(myslabP, N, steady=True, fmt='csc')
kP1 = eigenproblem(myPN, 'kappa', myslabP, nev=nev)
kP1.solve(algo=algo)
print('kP1={:5f}'.format(kP1.solution.eigvals[0]))

# S2
N = N+1
bc = 'Mark'
myslabS = Slab(M, xlayers, mats, [bc], G, N, 'FV')
mySN = NTE.SN(myslabS, N, steady=True, fmt='csc', BC=True)
kS2 = eigenproblem(mySN, 'kappa', myslabS, nev=nev)
kS2.solve(algo=algo)
print('kS2={:5f}'.format(kS2.solution.eigvals[0]))

phiD_g1 = kD.solution.get(group=1, moment=0, mode=0, normalise='peaktotalflux')
phiP1_g1 = kP1.solution.get(group=1, moment=0, mode=0, normalise='peaktotalflux')
phiS2_g1 = kS2.solution.get(group=1, moment=0, mode=0, normalise='peaktotalflux')


plt.plot(myslabD.mesh, phiD_g1, label='Diffusion')
plt.plot(myslabP.mesh, phiP1_g1, label='P1', linestyle=':')
plt.plot(myslabS.mesh, phiS2_g1, label='S2', linestyle='--')

if G > 1:
    phiD_g2 = kD.solution.get(group=2, moment=0, mode=0)
    phiP1_g2 = kP1.solution.get(group=2, moment=0, mode=0)
    phiS2_g2 = kS2.solution.get(group=2, moment=0, mode=0)
    
    plt.plot(myslabD.mesh, phiD_g2, label='Diffusion')
    plt.plot(myslabP.mesh, phiP1_g2, label='P1')
    plt.plot(myslabS.mesh, phiS2_g2, label='S2', linestyle='--')

# plt.plot(myslabSFD.mesh, phiSN1_FD, label='SN-FD', linestyle=':')
plt.axvline(x=0, c='k')
plt.legend()
plt.show()

# # get operators
# Lp1=myPN.L.todense()
# Sp1=myPN.S.todense()
# Fp1=myPN.F.todense()
# Rp1=myPN.R.todense()
# S0p1=myPN.S0.todense()
# F0p1=myPN.F0.todense()

# Ls2=mySN.L.todense()
# Ss2=mySN.S.todense()
# Fs2=mySN.F.todense()
# Rs2=mySN.R.todense()
# S0s2=mySN.S0.todense()
# F0s2=mySN.F0.todense()
