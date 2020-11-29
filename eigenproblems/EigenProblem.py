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
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


class eigenproblem():

    def __init__(self, nte, which):

        self.nS = nte.nS
        self.nE = nte.nE
        self.nA = nte.nA
        self.BC = nte.BC
        self.problem = which

    def _petsc(L, nev, what, P=None, which='LM', verbosity=False, sigma=0):

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

            # PC = PETSc.PC()  # PCFactorSetShiftType('NONZERO')
            # PC.setFactorSolverType('matsolverumfpack')
            # PC.setFactorShift(shift_type='nonzero', amount=0)
            # P.reorderForNonzeroDiagonal(atol=1E-12)

        if what in ['alpha', 'delta']:

            if which in ['SM', 'SR']:
                # E settings
                E = SLEPc.EPS().create()

                if P is not None:
                    E.setOperators(P, L)
                    E.setDimensions(nev=nev)
                    E.setProblemType(SLEPc.EPS.ProblemType.GNHEP)

                else:
                    E.setOperators(L)
                    E.setDimensions(nev=nev)
                    E.setProblemType(SLEPc.EPS.ProblemType.NHEP)

                E.setWhichEigenpairs(E.Which.TARGET_MAGNITUDE)
                if sigma is not None:
                    E.setTarget(sigma)

                st = E.getST()
                st.setType('sinvert')

            else:

                # E settings
                E = SLEPc.EPS().create()
                if P is not None:
                    E.setOperators(P, L)
                    E.setDimensions(nev=nev)
                    E.setProblemType(SLEPc.EPS.ProblemType.GNHEP)

                else:
                    E.setOperators(L)
                    E.setDimensions(nev=nev)
                    E.setProblemType(SLEPc.EPS.ProblemType.NHEP)

                if which == 'LM':
                    E.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_MAGNITUDE)

                elif which == 'LR':
                    E.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_REAL)

        elif what in ['gamma', 'kappa']:
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
        print('under develop')

    def plot(self, geom, moment, group, mode, ax=None, title=None, imag=False):
        yr, yi = eigenproblem.get(self, geom, moment, group, mode)
        x = geom.mesh if len(yr) == geom.NT else geom.stag_mesh
        ax = ax or plt.gca()
        if imag is False:
            plt.plot(x, yr)
        else:
            plt.plot(x, yi)
        ax.locator_params(nbins=8)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))

    def plotspectrum(self, loglog=False, ax=None, gaussplane=True):
        ax = ax or plt.gca()
        if gaussplane is True:
            plt.scatter(self.eigvals[1:].real, self.eigvals[1:].imag,
                        marker='o', color='red')
            # plot fundamental
            val, vect = eigenproblem.getfundamental(self)
            plt.scatter(val.real, val.imag, marker='*',
                        s=100, color='blue')
        else:
            plt.scatter(np.arange(0, len(self.eigvals.real)-1),
                        self.eigvals[1:].real, marker='o', color='red')
            # plot fundamental
            val, vect = eigenproblem.getfundamental(self)
            plt.scatter(0, val, marker='*', s=100, color='blue')

        if loglog is True:
            plt.yscale('symlog')
            plt.xscale('symlog')

        if self.problem == 'alphaprompt':
            label = 'alpha_p'
        elif self.problem == 'alphaprompt':
            label = 'alpha_d'
        else:
            label = self.problem

        plt.xlabel('$Re(\%s)$' % label)
        plt.ylabel('$Im(\%s)$' % label)
        # ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        # ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))

    def polarspectrum(self, ax=None):
        ax = ax or plt.gca()
        plt.polar(np.angle(self.eigvals[1:]), abs(self.eigvals[1:]),
                  marker='o', color='red')
        # plot fundamental
        val, vect = eigenproblem.getfundamental(self)
        plt.polar(np.angle(val), abs(val), marker='*', color='blue')

    def getfundamental(self):
        if self.problem in ['kappa', 'gamma']:
            idx = 0
        elif self.problem  == 'delta':
            # select real eigenvalues
            reals = self.eigvals[self.eigvals.imag == 0]
            # FIXME patch to find delta
            reals = reals[reals > 0]
            reals = reals[reals < 10]
            minreal = np.where(reals == reals.max())
            idx = np.where(self.eigvals == reals[minreal])[0]
        else:
            # select real eigenvalues
            reals = self.eigvals[self.eigvals.imag == 0]
            reals_abs = abs(self.eigvals[self.eigvals.imag == 0])
            minreal = np.where(reals_abs == reals_abs.min())
            idx = np.where(self.eigvals == reals[minreal])[0]

        return self.eigvals[idx], self.eigvect[:, idx]

    def get(self, geom, moment, group, mode):
        """
        Get spatial flux distribution for group, moment and spatial mode.

        Parameters
        ----------
        geom : object
            Geometry object.
        moment : int
            Moment number.
        group : int
            Group number.
        mode : int
            Eigen-mode number.

        Returns
        -------
        yr : ndarray
            Real part of the flux mode.
        yi : ndarray
            Imaginary part of the flux mode.

        """
        Neven = sum(np.arange(0, moment) % 2 == 0)
        Nodd = moment-Neven
        NT = geom.NT
        nE = self.nE
        iseven = 0 if moment % 2 else 1
        if Nodd >= 0:
            iS = Neven*(nE*NT)+Nodd*(nE*(NT-1))+(group-1)*((NT-1)+iseven)
        else:
            iS = group*(Neven*nE*NT)
        iE = iS+(NT-1)+iseven

        if mode == 0:
            _, vect = eigenproblem.getfundamental(self)
        else:
            vect = self.eigvect[:, mode]

        yr = np.real(vect[iS:iE])
        yi = np.imag(vect[iS:iE])
        return yr, yi
