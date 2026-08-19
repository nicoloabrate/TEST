"""
Author: N. Abrate.

File: test_ANL_6A1_reflected.py

Description: Numerical benchmark for a reflected 2G slab.
"""
import sys
sys.path.append('../../')
import pytest
from TEST.geometry import Slab
from TEST.models.NeutronTransportEquation import NTE
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
    myDiff = NTE(myslabD, "Diffusion", N=0, steady=True, fmt='csc')
    kD = eigenproblem(nte=myDiff, which='kappa', ge=myslabD, nev=nev)
    print('Elapsed time to setup Diffusion FD: {} s'.format(t.time()-t0))

    t0 = t.time()
    kD.solve(algo=algo, normalisation='peaktotalflux')
    print(f'Elapsed time, diffusion: {t.time()-t0:.3f} s')

    # display solution
    fig, ax = plt.subplots()
    myslabD.displaygeom()
    kD.solution.xplot(1)
    kD.solution.xplot(2)
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
    myPN = NTE(myslabP, "PN", N=N, steady=True, fmt='csc')
    kP1 = eigenproblem(nte=myPN, which='kappa', ge=myslabP, nev=1)
    print(f'Elapsed time to setup P{N} FD: {t.time()-t0:.3f}  s')

    t0 = t.time()
    kP1.solve(algo=algo, normalisation='peaktotalflux')
    print(f'Elapsed time, P{N}: {t.time()-t0:.3f} s')
    assert abs(1E5*(kP1.solution.eigvals[0]-kref)) < 800

    fig, ax = plt.subplots()
    myslabP.displaygeom()
    kP1.solution.xplot(1)
    kP1.solution.xplot(2)
    plt.title(f'P{N}')

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
    mySNFD = NTE(myslabSFD, "SN", N=N, steady=True, fmt='csc', BC=True)
    kS2FD = eigenproblem(nte=mySNFD, which='kappa', ge=myslabSFD, nev=1)
    print(f'Elapsed time to setup S{N} FD: {t.time()-t0} s')

    t0 = t.time()
    kS2FD.solve(algo=algo, normalisation='peaktotalflux')
    print(f'Elapsed time, S{N} FD: {t.time()-t0} s')
    assert abs(1E5*(kS2FD.solution.eigvals[0]-kref)) < 800

    fig, ax = plt.subplots()
    myslabSFD.displaygeom()
    kS2FD.solution.xplot(1)
    kS2FD.solution.xplot(2)
    plt.title(f'S{N} FD')

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
    mySNFV = NTE(myslabSFV, "SN", N=N, steady=True, fmt='csc', BC=True)
    kS2FV = eigenproblem(nte=mySNFV, which='kappa', ge=myslabSFV, nev=1)
    print(f'Elapsed time to setup S{N} FV: {t.time()-t0} s')

    t0 = t.time()
    kS2FV.solve(algo=algo, normalisation='peaktotalflux')
    print(f'Elapsed time, S{N} FV: {t.time()-t0} s')

    fig, ax = plt.subplots()
    myslabSFV.displaygeom()
    kS2FV.solution.xplot(1)
    kS2FV.solution.xplot(2)
    plt.title(f'S{N} FV')

    assert abs(1E5*(kS2FV.solution.eigvals[0]-kref)) < 800



if __name__ == '__main__':

    algo='eigs'
    M=[-21, -81, -21]
    xlayers = [-120, -80, 80, 120]
    mats=['ANL-6A1-R1-R3', 'ANL-6A1-R2', 'ANL-6A1-R1-R3']
    test_ANL_6A1_reflected_S2_FV(algo, M, xlayers, mats)
