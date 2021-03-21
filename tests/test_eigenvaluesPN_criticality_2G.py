"""
Author: N. Abrate.

File: test_eigenvaluesPN_1G.py

Description: Benchmark eigenvalues computed with PN approximation for a 
             critical system.
"""
import sys
sys.path.append('../../')
import time as t
from TEST.geometry import Slab
import TEST.NeutronTransportEquation as NTE
from TEST.eigenproblems.time import alpha
from TEST.eigenproblems.criticality import kappa
from TEST.eigenproblems.collision import gamma
from TEST.eigenproblems.streaming import delta


def test_PNcriticality(N):
    nev = 1
    M = 100
    G = 2
    bc = 'Mark'
    H = 1.795602
    
    matname = ['Pu239a']
    xlayers = [-H, H]
    # define geometry and mesh
    slab = Slab(M, xlayers, matname, [bc], G, N, 'FD')
    PN = NTE.PN(slab, N, steady=False, fmt='csc', prompt=True)
    
    a = alpha(slab, PN, nev=nev, verbosity=True)
    g = gamma(slab, PN, nev=nev, verbosity=True)
    d = delta(slab, PN, nev=nev, verbosity=True)
    k = kappa(slab, PN, nev=nev, verbosity=True)
    
    flxk1, _ = k.get(slab, 1, angle=0, mode=0)
    flxk2, _ = k.get(slab, 2, angle=0, mode=0)
    
    assert abs(k.eigvals[0]-1)<10 and abs(g.eigvals[0]-1)<10 and abs(d.eigvals[0]-1)<10

