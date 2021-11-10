"""
Author: N. Abrate.

File: test_ANL_6A1_reflected.py

Description: Numerical benchmark for a reflected 2G slab.
"""
import sys
sys.path.append('../../')
import pytest
from TEST.geometry import Slab
import TEST.models.NeutronTransportEquation as NTE
from TEST.models.EigenProblem import eigenproblem
import matplotlib.pyplot as plt
import time as t

@pytest.mark.parametrize("algo",['eigs', 'SLEPc'])
@pytest.mark.parametrize("M, xlayers, mats",
                          [([-21, -81, -21], [-120, -80, 80, 120],
                           ['ANL-6A1-R1-R3', 'ANL-6A1-R2', 'ANL-6A1-R1-R3'])])

def test_ANL_6A1_reflected_Diffusion(algo, M, xlayers, mats):

    kref = 0.9015507

    G = 2
    nev = 1

    # Diffusion
    N = 0
    bc = 'zero'
    myslabD = Slab(M, xlayers, mats, bc, G, N, 'FD')

    t0 = t.time()
    myDiff = NTE.Diffusion(myslabD, steady=True, fmt='csc')
    kD = eigenproblem(myDiff, 'kappa', myslabD, nev=nev)
    print('Elapsed time to setup PN FD: {} s'.format(t.time()-t0))

    t0 = t.time()
    kD.solve(algo=algo, normalisation='peaktotalflux')
    print('Elapsed time, diffusion: {} s'.format(t.time()-t0))

    # display solution
    fig, ax = plt.subplots()
    myslabD.displaygeom()
    kD.solution.plot(1)
    kD.solution.plot(2)
    plt.title('Diffusion')

    assert abs(1E5*(kD.solution.eigvals[0]-kref)) < 10


@pytest.mark.parametrize("algo",['eigs', 'SLEPc'])
@pytest.mark.parametrize("M, xlayers, mats",
                          [([-21, -81, -21], [-120, -80, 80, 120],
                           ['ANL-6A1-R1-R3', 'ANL-6A1-R2', 'ANL-6A1-R1-R3'])])

def test_ANL_6A1_reflected_P1(algo, M, xlayers, mats):

    kref = 0.9015507

    N = 1
    G = 2
    bc = 'Mark'
    myslabP = Slab(M, xlayers, mats, bc, G, N, 'FD')

    t0 = t.time()
    myPN = NTE.PN(myslabP, N, steady=True, fmt='csc')
    kP1 = eigenproblem(myPN, 'kappa', myslabP, nev=1)
    print('Elapsed time to setup PN FD: {} s'.format(t.time()-t0))

    t0 = t.time()
    kP1.solve(algo=algo, normalisation='peaktotalflux')
    print('Elapsed time, PN: {} s'.format(t.time()-t0))
    assert abs(1E5*(kP1.solution.eigvals[0]-kref)) < 800

    fig, ax = plt.subplots()
    myslabP.displaygeom()
    kP1.solution.plot(1)
    kP1.solution.plot(2)
    plt.title('P1')

@pytest.mark.parametrize("algo",['eigs', 'SLEPc'])
@pytest.mark.parametrize("M, xlayers, mats",
                          [([-21, -81, -21], [-120, -80, 80, 120],
                           ['ANL-6A1-R1-R3', 'ANL-6A1-R2', 'ANL-6A1-R1-R3'])])

def test_ANL_6A1_reflected_S2_FD(algo, M, xlayers, mats):

    kref = 0.9015507

    N = 2
    G = 2
    bc = 'Mark'
    myslabSFD = Slab(M, xlayers, mats, bc, G, N, 'FD')

    t0 = t.time()
    mySNFD = NTE.SN(myslabSFD, N, steady=True, fmt='csc', BC=True)
    kS2FD = eigenproblem(mySNFD, 'kappa', myslabSFD, nev=1)
    print('Elapsed time to setup SN FD: {} s'.format(t.time()-t0))

    t0 = t.time()
    kS2FD.solve(algo=algo, normalisation='peaktotalflux')
    print('Elapsed time, SN FD: {} s'.format(t.time()-t0))
    assert abs(1E5*(kS2FD.solution.eigvals[0]-kref)) < 800

    fig, ax = plt.subplots()
    myslabSFD.displaygeom()
    kS2FD.solution.plot(1)
    kS2FD.solution.plot(2)
    plt.title('S2 FD')

@pytest.mark.parametrize("algo",['eigs', 'SLEPc'])
@pytest.mark.parametrize("M, xlayers, mats",
                          [([-21, -81, -21], [-120, -80, 80, 120],
                           ['ANL-6A1-R1-R3', 'ANL-6A1-R2', 'ANL-6A1-R1-R3'])])

def test_ANL_6A1_reflected_S2_FV(algo, M, xlayers, mats):

    kref = 0.9015507

    N = 2
    G = 2
    bc = 'Mark'
    myslabSFV = Slab(M, xlayers, mats, bc, G, N, 'FV')

    t0 = t.time()
    mySNFV = NTE.SN(myslabSFV, N, steady=True, fmt='csc', BC=True)
    kS2FV = eigenproblem(mySNFV, 'kappa', myslabSFV, nev=1)
    print('Elapsed time to setup SN FV: {} s'.format(t.time()-t0))

    t0 = t.time()
    kS2FV.solve(algo=algo, normalisation='peaktotalflux')
    print('Elapsed time, SN FV: {} s'.format(t.time()-t0))
    print('kS2 FV={:5f}'.format(kS2FV.solution.eigvals[0]))

    fig, ax = plt.subplots()
    myslabSFV.displaygeom()
    kS2FV.solution.plot(1)
    kS2FV.solution.plot(2)
    plt.title('S2 FV')

    assert abs(1E5*(kS2FV.solution.eigvals[0]-kref)) < 800
