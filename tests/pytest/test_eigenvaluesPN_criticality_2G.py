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
from TEST.models.NeutronTransportEquation import NTE
from TEST.models.EigenProblem import eigenproblem


def test_PNcriticality():
    nev = 1
    M = 100
    G = 2
    N = 51
    bc = 'Mark'
    H = 1.795602

    matname = ['Pu239a']
    xlayers = [-H, H]
    # define geometry and mesh
    myslab = Slab(M, xlayers, matname, [bc], G, N, 'FD')
    myPN = NTE(myslab, "PN", N=N, steady=False, fmt='csc', prompt=True)

    a = eigenproblem(nte=myPN, which='alpha', ge=myslab, nev=nev)
    g = eigenproblem(nte=myPN, which='gamma', ge=myslab, nev=nev)
    g.solve(algo='SLEPc')
    d = eigenproblem(nte=myPN, which='delta', ge=myslab, nev=nev+5)
    d.solve(algo='SLEPc')
    k = eigenproblem(nte=myPN, which='kappa', ge=myslab, nev=nev)
    k.solve(algo='SLEPc')

    # flxk1, _ = k.solution.get(1, angle=0, mode=0)
    # flxk2, _ = k.solution.get(2, angle=0, mode=0)

    assert abs(k.solution.eigvals[0]-1)<10 and abs(g.solution.eigvals[0]-1)<10 and abs(d.solution.eigvals[0]-1)<10

def test_SNcriticality():
    nev = 1
    M = 100
    G = 2
    N = 51
    bc = 'Mark'
    H = 1.795602

    matname = ['Pu239a']
    xlayers = [-H, H]
    # define geometry and mesh
    myslab = Slab(M, xlayers, matname, [bc], G, N+1, 'FD')
    mySN = NTE(myslab, "SN", N=N+1, steady=False, fmt='csc', prompt=True)

    a = eigenproblem(nte=mySN, which='alpha', ge=myslab, nev=nev)
    g = eigenproblem(nte=mySN, which='gamma', ge=myslab, nev=nev)
    g.solve(algo='SLEPc')
    d = eigenproblem(nte=mySN, which='delta', ge=myslab, nev=nev+2)
    d.solve(algo='SLEPc')
    k = eigenproblem(nte=mySN, which='kappa', ge=myslab, nev=nev)
    k.solve(algo='SLEPc')

    # flxk1, _ = k.solution.get(1, angle=0, mode=0)
    # flxk2, _ = k.solution.get(2, angle=0, mode=0)

    assert abs(k.solution.eigvals[0]-1)<10 and abs(g.solution.eigvals[0]-1)<10 and abs(d.solution.eigvals[0]-1)<10
