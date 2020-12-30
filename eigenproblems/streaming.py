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

    def __init__(self, geom, nte, nev=1, algo='PETSc', verbosity=False,
                 normalization=None):

        super(delta, self).__init__(nte, 'delta')

        # define eigenproblem operators
        L = nte.L  # leakage operator
        B = nte.F+nte.S-nte.R  # material operator

        if algo == 'PETSc':

            try:
                start = t.time()
                eigvals, eigvect = eigenproblem._petsc(B, nev, 'delta',  P=L,
                                                       which='SR', sigma=1)
                eigvals = 1/eigvals
                end = t.time()

            except NameError:
                print('PETSc/SLEPc packages not available.')

                if L.format != 'csc':
                    L = L.tocsc()

                if B.format != 'csc':
                    B = B.tocsc()

                start = t.time()
                eigvals, eigvect = eigs(B, M=L, k=nev, which='LM')
                end = t.time()
                algo = 'eigs'

        elif algo == 'eigs':

            if L.format != 'csc':
                L = L.tocsc()

            if B.format != 'csc':
                B = B.tocsc()

            start = t.time()
            eigvals, eigvect = eigs(B, M=L, k=nev, which='LR', sigma=1)
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

        self.eigvals = eigvals
        self.eigvect = eigvect
