#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 10:17:07 2021

@author: nabrate
"""

import sys
sys.path.append('C:\\Users\\39346\\Documenti\\mycodes')
sys.path.append('/opt/programs_nabrate/mycodes/TEST_GTE')
import numpy as np
from TEST.geometry import Slab
import TEST.models.NeutronTransportEquation as NTE
from TEST.models.EigenProblem import eigenproblem
from TEST.models.GeneralizedEigenTheory import GET
import matplotlib.pyplot as plt

nev = 4
M = [-4, -5, -4]
G = 1
N = 0
# bc = 'Mark'
bc = 'zero'

H = 80
a = 40
matname = ['Fissile', 'Absorber', 'Fissile']
xlayers = [-H, -a, a, H]
# define geometry and mesh
myslab = Slab(M, xlayers, matname, [bc], G, N, 'FD')
# myPN = NTE.PN(myslab, N, steady=True, fmt='csc')
myPN = NTE.Diffusion(myslab, steady=True, fmt='csc')
# --- kappa eigenvalue
k = eigenproblem(myPN, 'kappa', myslab, nev=nev)
k.solve(algo='eig')
keff, _ = k.solution.getfundamental()
# plot
fig, ax = plt.subplots()
myslab.displaygeom()
k.solution.plot(1, mode=0, title='k: {:.5f}'.format((keff)))
plt.show()


Lp1=np.asarray(myPN.L.todense())
LInf=np.asarray(myPN.Linf.todense())
Sp1=np.asarray(myPN.S.todense())
Fp1=np.asarray(myPN.F.todense())
Rp1=np.asarray(myPN.R.todense())
Cp1=np.asarray(myPN.C.todense())
S0p1=np.asarray(myPN.S0.todense())
F0p1=np.asarray(myPN.F0.todense())


# analytic
A1=myslab.regions['Fissile'].Abs
D1=myslab.regions['Fissile'].Diffcoef
dx = myslab.dx
c1 = 2*D1/dx[0]**2