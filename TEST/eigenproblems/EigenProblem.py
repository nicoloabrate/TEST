"""
Author: N. Abrate.

File: EigenProblem.py

Description: Class that defines and solve the different eigenvalue problems
             defined in linear neutron transport theory.
"""
try:
    from petsc4py import PETSc
    from slepc4py import SLEPc

except ImportError:
    print('WARNING: PETSc/SLEPc packages not available. Computational speed' +
          ' may be seriously affected.')

import time as t
import scipy.sparse as scisparse
import numpy as np
import scipy.linalg as scilinalg
from numpy.linalg import norm


class eigenproblem():

    def __init__(self, nte, which):

        self.nS = nte.nS
        self.nE = nte.nE
        self.nA = nte.nA
        self.BC = nte.BC
        self.problem = which

    def petsc(L, nev, what, P=None, which='LM', verbosity=False):

        # BUG: set to 0 explicitly diagonal terms that does not appear (and thus are null)
        if L.format != 'csr':
            L = L.tocsr()

        # PETSc requires full diagonal
        diagL = L.diagonal()
        idL = np.array(np.where([diagL == 0])[1])
        # explicitly force 0 on diagonal
        L[idL, idL] = 0  # np.zeros((len(idL), 0))

        if P is not None:
            P = P.tocsr()
            diagP = P.diagonal()
            idP = np.array(np.where([diagP == 0])[1])
            P[idP, idP] = 0

        rows, cols = L.shape
        # create PETSc matrices
        L = PETSc.Mat().createAIJ(size=L.shape,
                                  csr=(L.indptr, L.indices, L.data))
        if P is not None:
            P = PETSc.Mat().createAIJ(size=P.shape,
                                      csr=(P.indptr, P.indices, P.data))

        if what == 'alpha':

            if which in ['SM', 'SR']:
                # create spectral transformation
                pc = PETSc.PC().create()
                pc.setType(pc.Type.BJACOBI)

                ksp = PETSc.KSP().create()
                ksp.setType(ksp.Type.PREONLY)
                ksp.setPC(pc)

                F = SLEPc.ST().create()
                F.setType(F.Type.SINVERT)
                F.setKSP(ksp)
                F.setShift(0.2)

                # E settings
                E = SLEPc.EPS().create()
                E.setST(F)

                if P is not None:
                    E.setOperators(P, L)
                    E.setDimensions(nev=nev)
                    E.setWhichEigenpairs(E.Which.LARGEST_MAGNITUDE)
                    E.setProblemType(SLEPc.EPS.ProblemType.GNHEP)

                else:
                    E.setOperators(L)
                    E.setDimensions(nev=nev)
                    E.setWhichEigenpairs(E.Which.LARGEST_MAGNITUDE)
                    E.setProblemType(SLEPc.EPS.ProblemType.NHEP)

            else:

                # E settings
                E = SLEPc.EPS().create()
                if P is not None:
                    E.setOperators(P, L)
                    E.setProblemType(SLEPc.EPS.ProblemType.GNHEP)

                else:
                    E.setOperators(L)
                    E.setProblemType(SLEPc.EPS.ProblemType.NHEP)

                if which == 'LM':
                    E.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_MAGNITUDE)

                elif which == 'LR':
                    E.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_REAL)

        elif what in ['gamma', 'delta', 'kappa']:
            E = SLEPc.EPS().create()
            E.setOperators(P, L)
            E.setWhichEigenpairs(E.Which.LARGEST_MAGNITUDE)
            E.setDimensions(nev=nev)

        E.setFromOptions()

        start = t.time()
        E.solve()
        end = t.time()

        if verbosity is not None:
            print("ELAPSED TIME: %f [s]" % (end-start))

        vr, vi = L.getVecs()

        vals = []
        vecs = []
        eigvect = np.full((rows, max(nev, E.getConverged())), np.nan)
        for iE in range(E.getConverged()):
            val = E.getEigenpair(iE, vr, vi)
            vals.append(val)
            vecs = [complex(vr0, wi0) for vr0, wi0 in zip(vr.getArray(),
                                                          vi.getArray())]
            eigvect[:, iE] = np.asarray(vecs, dtype=np.complex).T

        eigvals = np.asarray(vals)
        return eigvals, eigvect

    def normalize(eigvect, which=None):
        """Normalize eigenvectors according to a user-defined criterion."""
        print('develop')
