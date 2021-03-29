"""
Author: N. Abrate.

File: test_eigenvaluesPN_1G.py

Description: Benchmark eigenvalues computed with PN approximation.
"""
import sys
sys.path.append('../../')
import pytest
from TEST.geometry import Slab
import TEST.NeutronTransportEquation as NTE
from TEST.eigenproblems.EigenProblem import eigenproblem


@pytest.mark.filterwarnings('ignore::DeprecationWarning:SparseEfficiencyWarning')
@pytest.mark.parametrize("algo",['eigs', 'PETSc'])
@pytest.mark.parametrize("H, N, ref",
                         [(8, 3, 4.225105), (8, 7, 4.229065), (8, 15, 4.229840),
                          (8, 31, 4.230020), (8, 255, 4.230078), (1, 3, 1.129075),
                          (1, 7, 1.209376), (1, 15, 1.223479), (1, 31, 1.225737),
                          (1, 255, 1.226406)])
def test_Modak_kappa0_1G(H, N, ref, algo):
    """
    Benchmark based on the eigenvalues taken from ``Modak, R. S., D. C. Sahni, 
    and S. D. Paranjape. 1995. "Evaluation of higher k-eigenvalues of the 
    neutron transport equation by SN method" Ann. Nucl. Energy 22 (6):359–66.``
    Test for one-group PN and SN modules.

    Returns
    -------
    None.

    """
    nev = 1
    M = 100
    G = 1
    bc = 'Mark'
    xlayers = [0, H]
    # define geometry and mesh
    myslab = Slab(M, xlayers, ['Modak'], [bc], G, N, 'FD')
    myPN = NTE.PN(myslab, N, steady=True, fmt='csc')
    k1 = eigenproblem(myPN, 'kappa', myslab, nev=nev)
    k1.solve(verbosity=True, algo=algo)
    assert abs(k1.solution.eigvals[0]-ref)*1E5 < 5

@pytest.mark.parametrize("algo",['eigs', 'PETSc'])
def test_Modak_kappa0_aniso_1G(algo):
    """
    Benchmark based on the eigenvalues taken from ``Modak, R. S., D. C. Sahni, 
    and S. D. Paranjape. 1995. "Evaluation of higher k-eigenvalues of the 
    neutron transport equation by SN method" Ann. Nucl. Energy 22 (6):359–66.``
    Test for one-group PN and SN modules.

    Returns
    -------
    None.

    """
    ref = 4.0653591
    nev = 1
    M = 100
    G = 1
    H = 8
    N = 63
    bc = 'Mark'
    xlayers = [0, H]
    # define geometry and mesh
    myslab = Slab(M, xlayers, ['ModakAni'], [bc], G, N, 'FD')
    myPN = NTE.PN(myslab, N, steady=True, fmt='csc')
    k1 = eigenproblem(myPN, 'kappa', myslab, nev=nev)
    k1.solve(verbosity=True, algo=algo)
    assert abs(k1.solution.eigvals[0]-ref)*1E5 < 5

@pytest.mark.parametrize("algo",['eigs', 'PETSc'])
@pytest.mark.filterwarnings('ignore::DeprecationWarning:SparseEfficiencyWarning')
def test_Modak_kappa_higher_1G(algo):
    """
    Benchmark based on the eigenvalues taken from ``Modak, R. S., D. C. Sahni, 
    and S. D. Paranjape. 1995. "Evaluation of higher k-eigenvalues of the 
    neutron transport equation by SN method" Ann. Nucl. Energy 22 (6):359–66.``
    Test for one-group PN and SN modules.

    Returns
    -------
    None.

    """
    ref = [4.2300510, 2.041626, 1.13835, 0.75610, 0.5583, 0.440]
    tol = [10, 30, 80, 90, 100, 150]
    nev = 10
    M = 100
    G = 1
    H = 8
    N = 63
    bc = 'Mark'
    xlayers = [0, H]
    # define geometry and mesh
    myslab = Slab(M, xlayers, ['Modak'], [bc], G, N, 'FD')
    myPN = NTE.PN(myslab, N, steady=True, fmt='csc')
    k1 = eigenproblem(myPN, 'kappa', myslab, nev=nev)
    k1.solve(verbosity=True, algo=algo)
    for i, k in enumerate(k1.solution.eigvals[0::2]):
        assert abs(k-ref[i])*1E5 < tol[i]

@pytest.mark.parametrize("algo",['eigs', 'PETSc'])
@pytest.mark.filterwarnings('ignore::DeprecationWarning:SparseEfficiencyWarning')
def test_Modak_gamma_higher_1G(algo):
    """
    Benchmark based on the eigenvalues taken from ``Modak, R. S., D. C. Sahni, 
    and S. D. Paranjape. 1995. "Evaluation of higher k-eigenvalues of the 
    neutron transport equation by SN method" Ann. Nucl. Energy 22 (6):359–66.``
    Test for one-group PN and SN modules.

    Returns
    -------
    None.

    """
    ref = [1.0364038, 1.289806, 1.678464, 2.12257]
    tol = [5, 10, 50, 150]
    nev = 7
    M = 100
    G = 1
    H = 8
    N = 63
    bc = 'Mark'
    xlayers = [0, H]
    # define geometry and mesh
    myslab = Slab(M, xlayers, ['Modak'], [bc], G, N, 'FD')
    myPN = NTE.PN(myslab, N, steady=True, fmt='csc')
    g1 = eigenproblem(myPN, 'gamma', myslab, nev=nev)
    g1.solve(verbosity=True, algo=algo)
    for i, g in enumerate(g1.solution.eigvals[0::2]):
        g = 1/g*(1.8)  # c=(XS_S-NU*XS_F)/XS_T/gamma
        assert abs(g-ref[i])*1E5 < tol[i]

@pytest.mark.parametrize("algo",['eigs', 'PETSc'])
@pytest.mark.parametrize("N, ref, tol",
                         [(7, [-2.53782E-02, -1.03353E-01, -2.38497E-01, -4.43506E-01],
                              [5, 50, 50, 100]),
                          (15, [-2.53493E-02, -1.03642E-01, -2.37487E-01, -4.51828E-01],
                               [5, 70, 70, 1500])])
def test_Modak_alpha_higher_1G(N, ref, tol, algo):
    """
    Benchmark based on the eigenvalues taken from ``Modak, R. S., and A. Gupta.
    2003. "A simple scheme for the direct evaluation of time-eigen- values of 
    neutron transport equation". Ann. Nucl. Energy 30 (2):211–22.``
    Test for one-group PN and SN modules.

    Returns
    -------
    None.

    """
    nev = 4
    M = 50
    G = 1
    H = 10
    bc = 'Mark'
    xlayers = [0, H]
    # define geometry and mesh
    myslab = Slab(M, xlayers, ['Dahl'], [bc], G, N, 'FD')
    myPN = NTE.PN(myslab, N, steady=False, prompt=True, fmt='csc')
    a1 = eigenproblem(myPN, 'alpha', myslab, nev=nev)
    a1.solve(verbosity=True, algo=algo)
    for i, a in enumerate(a1.solution.eigvals):
        assert abs(a-ref[i])*1E5 < tol[i]