"""
Author: N. Abrate.

File: gamma.py

Description: collision eigenvalue problem class.
"""
import time as t
from scipy.sparse.linalg import eigs
import numpy as np
from scipy.linalg import eig
from .EigenProblem import eigenproblem


class gamma(eigenproblem):

    def __init__(self, geom, nte, nev=1, algo='PETSc', verbosity=False,
                 normalization=None, tol=1E-8):

        super(gamma, self).__init__(nte, 'gamma')
        res = None
        # define kappa eigenproblem operators
        if nev == 0:  # kappa infinite
            L = nte.R  # no leakage, infinite medium
            nev = 1
        else:
            L = nte.L + nte.R  # destruction operator

        P = nte.F + nte.S  # multiplication operator

        if algo == 'PETSc':

            try:
                start = t.time()
                eigvals, eigvect = eigenproblem._petsc(L, nev, 'gamma', P=P,
                                                      which='LM', tol=1E-8)
                end = t.time()

            except NameError:
                print('PETSc/SLEPc packages not available.')

                if L.format != 'csc':
                    L = L.tocsc()

                if P.format != 'csc':
                    P = P.tocsc()

                start = t.time()
                eigvals, eigvect = eigs(P, M=L, k=nev, which='LM')
                end = t.time()
                algo = 'eigs'

        elif algo == 'eigs':

            if L.format != 'csc':
                L = L.tocsc()

            if P.format != 'csc':
                P = P.tocsc()

            start = t.time()
            eigvals, eigvect = eigs(P, M=L, k=nev, which='LM')
            end = t.time()

        elif algo == 'eig':
            start = t.time()
            eigvals, eigvect = eig(P.todense(), L.todense())
            end = t.time()

        else:
            raise OSError('%s algorithm is unavailable!' % algo)

        if verbosity is True and algo != 'PETSc':
            print("ELAPSED TIME: %f [s]" % (end-start))

        self.algo = algo

        # sort eigenvalues
        idx = eigvals.argsort()[::-1]
        eigvals = eigvals[idx]
        eigvect = eigvect[:, idx]

        # force sign consistency
        signs = np.sign(eigvect[1, :])  # sign of 2nd row to avoid BCs
        eigvect = np.conj(signs)*eigvect

        # normalize eigenvectors
        for iv, v in enumerate(eigvect.T):
            eigvect[:, iv] = v/np.linalg.norm(v)

        # FIXME call balance function
        if res is not None:
             self.residual = res
        self.eigvals = eigvals
        self.eigvect = eigvect
