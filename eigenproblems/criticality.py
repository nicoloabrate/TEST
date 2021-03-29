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

    def __init__(self, geom, nte, nev=1):

        super(kappa, self).__init__(nte, 'kappa')
        # define kappa eigenproblem operators
        if nev == 0 or self.BC is False:  # kappa infinite
            L = nte.R-nte.S  # no leakage, infinite medium
            nev = 1
        else:
            L = nte.L+nte.R-nte.S  # destruction operator

        F = nte.F  # multiplication operator

        self.nev = nev
        if 2*nev+1 >= L.shape[0]:
            raise OSError('Too many eigenvalues required! 2*nev+1 should be < operator rank')

        self.A = L
        self.B = F

    def solve(self, algo='PETSc', verbosity=False, normalization=None,
              tol=1E-14, monitor=False, phasespace=None):

        L = self.A
        F = self.B
        res = None
        if algo == 'PETSc':

            try:
                start = t.time()
                eigvals, eigvect, res = eigenproblem._petsc(L, self.nev, 'kappa',  P=F,
                                                            which='LM', tol=tol,
                                                            monitor=monitor)
                end = t.time()

            except NameError:
                print('PETSc/SLEPc packages not available.')

                if L.format != 'csc':
                    L = L.tocsc()

                if F.format != 'csc':
                    F = F.tocsc()

                start = t.time()
                eigvals, eigvect = eigs(F, M=L, k=self.nev, which='LM')
                end = t.time()
                algo = 'eigs'

        elif algo == 'eigs':

            if L.format != 'csc':
                L = L.tocsc()

            if F.format != 'csc':
                F = F.tocsc()

            start = t.time()
            eigvals, eigvect = eigs(F, M=L, k=self.nev, which='LM')
            end = t.time()

        elif algo == 'eig':
            start = t.time()
            eigvals, eigvect = eig(F.todense(), L.todense())
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

        self.eigvals = eigvals[0:self.nev]
        # convert to np.float64 if imaginary part is null
        if np.iscomplex(eigvect[:, 0:self.nev]).sum() == 0:
            self.eigvect = eigvect[:, 0:self.nev].real
        else:
            self.eigvect = eigvect[:, 0:self.nev]

        # normalize eigenvectors
        self.normalize(phasespace=phasespace)
        # FIXME call balance function

        if res is not None:
             self.residual = res

