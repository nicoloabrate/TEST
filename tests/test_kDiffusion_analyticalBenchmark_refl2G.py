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
import TEST.NeutronTransportEquation as NTE
import TEST.AdjointTransportEquation as ATE
from TEST.eigenproblems.EigenProblem import eigenproblem

@pytest.mark.filterwarnings('ignore::DeprecationWarning:SparseEfficiencyWarning')
@pytest.mark.parametrize("algo",['eigs', 'PETSc'])
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
    bc = 'zero'
    xlayers = [-R, -H, H, R]
    # define geometry and mesh
    myslab = Slab(M, xlayers, [matrefl, 'MontagniniFuel', matrefl], [bc], G, N, 'FD')
    myPN = NTE.Diffusion(myslab, steady=True, fmt='csc')
    k1 = eigenproblem(myPN, 'kappa', myslab, nev=nev)
    k1.solve()
    assert abs(k1.solution.eigvals[0]-ref)*1E5 < 1


@pytest.mark.filterwarnings('ignore::DeprecationWarning:SparseEfficiencyWarning')
@pytest.mark.parametrize("algo",['eigs', 'PETSc'])
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
    bc = 'zero'
    xlayers = [-R, -H, H, R]
    # define geometry and mesh
    myslab = Slab(M, xlayers, [matrefl, 'MontagniniFuel', matrefl], [bc], G, N, 'FD')
    myPN = ATE.Diffusion(myslab, steady=True, fmt='csc')
    k1 = eigenproblem(myPN, 'kappa', myslab, nev=nev)
    k1.solve()
    assert abs(k1.solution.eigvals[0]-ref)*1E5 < 1
