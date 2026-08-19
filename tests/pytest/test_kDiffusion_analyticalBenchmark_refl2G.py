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
from TEST.models.NeutronTransportEquation import NTE
from TEST.models.EigenProblem import eigenproblem

@pytest.mark.filterwarnings('ignore::DeprecationWarning:SparseEfficiencyWarning')
@pytest.mark.parametrize("algo",['eigs', 'SLEPc'])
@pytest.mark.parametrize("H, R, G, matrefl, ref", [(30, 70, 2, 'MontagniniReflector3', 1.004241348107076),
                                                   (40, 100, 2, 'MontagniniReflector2', 1.045766651960480),
                                                   (40, 100, 2, 'MontagniniReflector3', 1.020903109926185)])
def test_Diffusion_kappa0(H, R, G, matrefl, ref, algo):
    """
    Analytical benchmark for a reflected 2G slab.

    Returns
    -------
    None.

    """
    nev = 1
    M = 50
    N = 0
    bc = 'Mark'
    xlayers = [-R, -H, H, R]
    # Diffusion
    myslab = Slab(M, xlayers, [matrefl, 'MontagniniFuel', matrefl], bc, G, N, 'FD')
    myD = NTE(myslab, "Diffusion", N=N, steady=True, fmt='csc')
    kD = eigenproblem(nte=myD, which='kappa', ge=myslab, nev=nev)
    kD.solve(algo=algo)
    dk_Diff = abs(kD.solution.eigvals[0] - ref)*1E5
    print(f"Diffusion: k = {kD.solution.eigvals[0]:.8f}, dk = {dk_Diff:.3f} pcm")
    assert dk_Diff < 5
    # P1
    N = 1
    bc = 'Mark'
    myslab_P1 = Slab(M, xlayers, [matrefl, 'MontagniniFuel', matrefl], bc, G, N, 'FD')
    myPN = NTE(myslab_P1, "PN", N=N, steady=True, fmt='csc')
    kP1 = eigenproblem(nte=myPN, which='kappa', ge=myslab_P1, nev=nev)
    kP1.solve(algo=algo)
    dk_P1 = abs(kP1.solution.eigvals[0].real - ref)*1E5
    print(f"P1: k = {kP1.solution.eigvals[0]:.8f}, dk = {dk_P1:.3f} pcm")
    assert dk_P1 < 5


@pytest.mark.filterwarnings('ignore::DeprecationWarning:SparseEfficiencyWarning')
@pytest.mark.parametrize("algo",['eigs', 'SLEPc'])
@pytest.mark.parametrize("H, R, G, matrefl, ref", [(30, 70, 2, 'MontagniniReflector3', 1.004241348107076),
                                                   (40, 100, 2, 'MontagniniReflector2', 1.045766651960480),
                                                   (40, 100, 2, 'MontagniniReflector3', 1.020903109926185)])
def test_DiffusionAdjoint_kappa0(H, R, G, matrefl, ref, algo):
    """
    Analytical benchmark for a reflected 2G slab.

    Returns
    -------
    None.

    """
    nev = 1
    M = 50
    N = 0
    bc = 'Mark'
    xlayers = [-R, -H, H, R]
    # define geometry and mesh
    myslab = Slab(M, xlayers, [matrefl, 'MontagniniFuel', matrefl], bc, G, N, 'FD')
    myPN = NTE(myslab, "Diffusion", N=N, steady=True, fmt='csc', adjoint=True)
    k1 = eigenproblem(nte=myPN, which='kappa', ge=myslab, nev=nev)
    k1.solve()
    assert abs(k1.solution.eigvals[0]-ref)*1E5 < 5

