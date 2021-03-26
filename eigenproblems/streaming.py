"""
Author: N. Abrate.

File: delta.py

Description: streaming eigenvalue problem class.
"""
import time as t
from scipy.sparse.linalg import eigs
import numpy as np
from scipy.linalg import eig
from .EigenProblem import eigenproblem

class delta(eigenproblem):

    def __init__(self, geom, nte, nev=1):

        super(delta, self).__init__(nte, 'delta')
        # define eigenproblem operators
        L = nte.L  # leakage operator
        B = nte.F+nte.S-nte.R  # material operator
        self.A = L
        self.B = B
        self.nev = nev

    def solve(self, algo='PETSc', verbosity=False, normalization=None,
              tol=1E-8, monitor=False):

        L = self.A
        B = self.B
        res = None
        if algo == 'PETSc':

            try:
                start = t.time()
                eigvals, eigvect, res = eigenproblem._petsc(B, self.nev, 'delta',  P=L,
                                                            which='SR', sigma=1,
                                                            tol=tol,
                                                            monitor=monitor)
                eigvals = 1/eigvals
                end = t.time()

            except NameError:
                print('PETSc/SLEPc packages not available.')

                if L.format != 'csc':
                    L = L.tocsc()

                if B.format != 'csc':
                    B = B.tocsc()

                start = t.time()
                eigvals, eigvect = eigs(B, M=L, k=self.nev, which='LM')
                end = t.time()
                algo = 'eigs'

        elif algo == 'eigs':

            if L.format != 'csc':
                L = L.tocsc()

            if B.format != 'csc':
                B = B.tocsc()

            start = t.time()
            eigvals, eigvect = eigs(B, M=L, k=self.nev, which='LR', sigma=1)
            end = t.time()

        elif algo == 'eig':
            start = t.time()
            eigvals, eigvect = eig(L.todense(), B.todense())
            eigvals = 1/eigvals
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

        # FIXME call balance function

        if res is not None:
             self.residual = res

        self.eigvals = eigvals
        self.eigvect = eigvect
