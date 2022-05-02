import os
import numpy as np
import shutil as sh
import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
from TEST.geometry import Slab
from TEST.material import Mix
from TEST.models.NeutronTransportEquation import NTE
from TEST.models.EigenProblem import eigenproblem
from TEST.models.GeneralizedEigenvalueTheory import GET
from TEST.utils import get_energy_grid

# general code settings
nev = 1
Nx = 10
G = '1G'
N = 1
algo = 'SLEPc'
model = f'P{N}'
BC = 'Marshak'
H = 1
# define phase space and core geometry
xlayers = [-H, H]
mat = 'Pu239'

# collapse data if needed
datapath = None
NG = 1

slab = Slab(-Nx, xlayers, mat, BC, G, N, 'FD', datapath=datapath)
# %%
mat=slab.regions['Pu239']
mat.Diffcoef = 1/3/mat.Tot
slab.perturb({"Nubar": {"howmuch": [1/1.1250616548935823-1], 'where': [(-H, H)]}})

#  check k eigenproblem response
myP = NTE(slab, model, steady=True)
k = eigenproblem(nte=myP, which='kappa', ge=slab, nev=4)
k.solve(algo=algo)
kpsN = k.solution
keff, _ = kpsN.getfundamental()
print(f'keff is {kpsN.eigvals[0].real:10f}')    

#  check k eigenproblem response
myD = NTE(slab, "Diffusion", steady=True)
kD = eigenproblem(nte=myD, which='kappa', ge=slab, nev=4)
kD.solve(algo=algo)
kpsD = kD.solution
keffD, _ = kpsD.getfundamental()
print(f'keff with Diffusion is {kpsD.eigvals[0]:10f}')    

algo = 'eig'

# Diffusion
myD = NTE(slab, "Diffusion", steady=0, prompt=0)
oD = eigenproblem(nte=myD, which='omega', ge=slab, nev=4, generalisedTime=False)
oD.solve(algo=algo)
opsD = oD.solution
omegaD, _ = opsD.getfundamental()
print(f'Omega with diffusion is {omegaD}')    

# P1
myP = NTE(slab, "P1", steady=0, prompt=0)
oP = eigenproblem(nte=myP, which='omega', ge=slab, nev=4)
oP.solve(algo=algo)
opsP = oP.solution
omegaP, _ = opsP.getfundamental()
print(f'Omega with P1 is {omegaP}')    
