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
import TEST.models.NeutronTransportEquation as NTE
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
    myPN = NTE.PN(myslab, N, steady=False, fmt='csc', prompt=True)

    a = eigenproblem(myPN, 'alpha', myslab, nev=nev)
    g = eigenproblem(myPN, 'gamma', myslab, nev=nev)
    g.solve(algo='SLEPc')
    d = eigenproblem(myPN, 'delta', myslab, nev=nev+5)
    d.solve(algo='SLEPc')
    k = eigenproblem(myPN, 'kappa', myslab, nev=nev)
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
    mySN = NTE.SN(myslab, N+1, steady=False, fmt='csc', prompt=True)

    a = eigenproblem(mySN, 'alpha', myslab, nev=nev)
    g = eigenproblem(mySN, 'gamma', myslab, nev=nev)
    g.solve(algo='SLEPc')
    d = eigenproblem(mySN, 'delta', myslab, nev=nev+2)
    d.solve(algo='SLEPc')
    k = eigenproblem(mySN, 'kappa', myslab, nev=nev)
    k.solve(algo='SLEPc')

    # flxk1, _ = k.solution.get(1, angle=0, mode=0)
    # flxk2, _ = k.solution.get(2, angle=0, mode=0)

    assert abs(k.solution.eigvals[0]-1)<10 and abs(g.solution.eigvals[0]-1)<10 and abs(d.solution.eigvals[0]-1)<10
