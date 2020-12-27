"""
Author: N. Abrate.

File: kappa.py

Description: kappa eigenvalue problem class.
"""
import time as t
from scipy.sparse.linalg import eigs
import numpy as np
from scipy.linalg import eig
from .EigenProblem import eigenproblem

class kappa(eigenproblem):

    def __init__(self, geom, nte, nev=1, algo='PETSc', verbosity=None,
                 normalization=None):

        super(kappa, self).__init__(nte, 'kappa')

        # define kappa eigenproblem operators
        if nev == 0:  # kappa infinite
            L = nte.R-nte.S  # no leakage, infinite medium
            nev = 1
        else:
            L = nte.L+nte.R-nte.S  # destruction operator

        F = nte.F  # multiplication operator

        if algo == 'PETSc':

            try:
                start = t.time()
                eigvals, eigvect = eigenproblem._petsc(L, nev, 'kappa',  P=F,
                                                       which='LM')
                end = t.time()

            except NameError:
                print('PETSc/SLEPc packages not available.')

                if L.format != 'csc':
                    L = L.tocsc()

                if F.format != 'csc':
                    F = F.tocsc()

                start = t.time()
                eigvals, eigvect = eigs(F, M=L, k=nev, which='LM')
                end = t.time()
                algo = 'eigs'

        elif algo == 'eigs':

            if L.format != 'csc':
                L = L.tocsc()

            if F.format != 'csc':
                F = F.tocsc()

            start = t.time()
            eigvals, eigvect = eigs(F, M=L, k=nev, which='LM')
            end = t.time()

        elif algo == 'eig':
            start = t.time()
            eigvals, eigvect = eig(F.todense(), L.todense())
            end = t.time()

        else:
            raise OSError('%s algorithm is unavailable!' % algo)

        if verbosity is not None and algo != 'PETSc':
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

        self.eigvals = eigvals
        self.eigvect = eigvect
