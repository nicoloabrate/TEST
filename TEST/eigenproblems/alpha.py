"""
Author: N. Abrate.

File: alpha.py

Description: alpha eigenvalue problem class.
"""
import time as t
from scipy.sparse.linalg import eigs, inv
import numpy as np
from scipy.linalg import eig
from .EigenProblem import eigenproblem


class alphaprompt(eigenproblem):

    def __init__(self, geom, nte, nev=1, algo='PETSc', verbosity=None,
                 normalization=None, which='SM', generalized=False):

        super(alphaprompt, self).__init__(nte, 'alphaprompt')
        # nev = min(nev + 10, nte.L.shape[0]-5)
        # compute maximum velocity and normalise equation
        vmax = -1/(np.min(-nte.T))

        # define alpha prompt eigenproblem operators
        if generalized is False:

            # invert time operator (velocity reciprocal)
            invT = inv(nte.T)/vmax

            if nev == 0:  # kappa infinite
                B = nte.S+nte.F-nte.R  # no leakage, infinite medium
                nev = 1
            else:
                B = nte.S+nte.F-nte.R-nte.L  # destruction operator

            B = np.dot(invT, B)
            T = None

        else:

            if nev == 0:  # kappa infinite
                B = nte.S+nte.F-nte.R  # no leakage, infinite medium
                nev = 1

            else:
                B = nte.S+nte.F-nte.R-nte.L  # destruction operator

            B = B/vmax
            T = nte.T/vmax  # multiplication operator

        if algo == 'PETSc':

            try:
                start = t.time()
                eigvals, eigvect = eigenproblem.petsc(B, nev, 'alpha', P=T,
                                                      which=which,
                                                      verbosity=True)
                end = t.time()

            except NameError:
                print('PETSc/SLEPc packages not available.')

                if B.format != 'csc':
                    B = B.tocsc()

                if T.format != 'csc':
                    T = T.tocsc()

                if which == 'SM':
                    which = 'LM'  # ARPACK much faster with shift-invert
                    sigma = 0
                else:
                    sigma = 0

                start = t.time()
                eigvals, eigvect = eigs(B, M=T, k=nev, which=which,
                                        sigma=sigma)
                end = t.time()
                algo = 'eigs'

        elif algo == 'eigs':

            if B.format != 'csc':
                B = B.tocsc()

            if T is not None:
                if T.format != 'csc':
                    T = T.tocsc()

            if which == 'SM':
                which = 'LM'  # ARPACK much faster with shift-invert
                sigma = 0
            else:
                sigma = 0

            start = t.time()
            eigvals, eigvect = eigs(B, M=T, k=nev, which=which, sigma=sigma)
            end = t.time()

        elif algo == 'eig':

            if T is None:
                start = t.time()
                eigvals, eigvect = eig(B.todense())
                end = t.time()
            else:
                start = t.time()
                eigvals, eigvect = eig(B.todense(), T.todense())
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
