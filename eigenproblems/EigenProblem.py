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
from copy import copy
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import ConnectionPatch


class eigenproblem():

    def __init__(self, nte, which):

        self.nS = nte.nS
        self.nE = nte.nE
        self.nA = nte.nA
        self.BC = nte.BC
        self.problem = which

    def _petsc(L, nev, what, P=None, which='LM', verbosity=False, sigma=0, tol=1E-8):

        start = t.time()
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

        end = t.time()

        if verbosity is True:
            print("ELAPSED TIME (PETSc setup): %f [s]" % (end-start))

        E.setTolerances(tol=tol) #
        start = t.time()
        E.setFromOptions()
        E.solve()
        end = t.time()

        if verbosity is True:
            print("ELAPSED TIME (PETSc solution): %f [s]" % (end-start))

        vr, vi = L.getVecs()

        vals = []
        vecs = []
        err = []
        eigvect = np.full((rows, max(nev, E.getConverged())), np.nan)
        for iE in range(E.getConverged()):
            val = E.getEigenpair(iE, vr, vi)
            vals.append(val)
            err.append(E.computeError(iE))
            vecs = [complex(vr0, wi0) for vr0, wi0 in zip(vr.getArray(),
                                                          vi.getArray())]
            eigvect[:, iE] = np.asarray(vecs, dtype=np.complex).T

        eigvals = np.asarray(vals)
        res = np.asarray(err)
        return eigvals, eigvect, res

    def issingular(A):
        A = A.todense()
        shapecheck = A.shape[0] == A.shape[1]
        rankcheck = np.linalg.matrix_rank(A) == A.shape[0]
        isinvertible = shapecheck and rankcheck
        return ~isinvertible

    def normalize(eigvect, which=None):
        """Normalize eigenvectors according to a user-defined criterion."""
        print('under develop')

    def plot(self, geom, group, moment=0, mode=0, family=0, precursors=False,
             ax=None, title=None, imag=False, **kwargs):

        yr, yi = eigenproblem.get(self, geom, group, moment=moment, mode=mode,
                                  family=family, precursors=precursors)
        x = geom.mesh if len(yr) == geom.NT else geom.stag_mesh
        ax = ax or plt.gca()

        if imag is False:
            plt.plot(x, yr, **kwargs)
        else:
            plt.plot(x, yi, **kwargs)

        ax.locator_params(nbins=8)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))

    def plotspectrum(self, loglog=False, gaussplane=True, geom=None,
                     grid=True, ylims=None, threshold=None, subplt=False):

        if self.problem == 'omega' and geom is not None:
            lambdas = geom.getxs('lambda')
            subplt = False if subplt is False else True
        else:
            lambdas = None

        if subplt is True:
            fig = plt.figure(figsize=(6.4*2, 4.8))
            sub1 = fig.add_subplot(1, 2, 1)
            sub2 = fig.add_subplot(1, 2, 2)
        else:
            fig = plt.figure()
            sub1 = fig.add_subplot(1, 1, 1)
            subplt = False

        val, vect = eigenproblem.getfundamental(self, lambdas)
        evals = np.delete(self.eigvals, np.where(self.eigvals == val))

        show = np.nan if threshold is not None and abs(val) > threshold else 1

        if threshold is not None:
            evals = evals[abs(evals) < threshold]

        if gaussplane is True:
            sub1.scatter(evals.real, evals.imag, marker='o', color='red')
            # plot fundamental
            sub1.scatter(show*val.real, show*val.imag, marker='*', s=100,
                         color='blue')
        else:
            sub1.scatter(np.arange(0, len(evals.real)-1), evals.real,
                         marker='o', color='red')
            # plot fundamental
            sub1.scatter(0, show*val, marker='*', s=100, color='blue')

        if self.problem == 'alpha':
            label = 'alpha'
        elif self.problem == 'omega':
            label = 'omega'
        else:
            label = self.problem

        if self.problem == 'alpha' or self.problem == 'omega':
            sub1.set_xlabel('$Re(\%s) ~[s^{-1}]$' % label)
            sub1.set_ylabel('$Im(\%s) ~[s^{-1}]$' % label)
        else:
            sub1.set_xlabel('$Re(\%s)$' % label)
            sub1.set_ylabel('$Im(\%s)$' % label)

        if ylims is None:
            sub1.set_ylim([min(self.eigvals.imag)*1.1, max(self.eigvals.imag)*1.1])
        else:
            sub1.set_ylim(ylims)

        if loglog is True:
            sub1.set_yscale('symlog')
            sub1.set_xscale('symlog')
        else:
            sub1.ticklabel_format(axis='x', scilimits=[-5, 5])
            sub1.ticklabel_format(axis='y', scilimits=[-5, 5])

        if grid is True:
            sub1.grid(alpha=0.2)

        if subplt is True:
            minl, maxl = min(-lambdas[:, 0]), max(-lambdas[:, 0])
            miny, maxy = min(self.eigvals.imag), max(self.eigvals.imag)
            minx, maxx = min(self.eigvals.real), max(self.eigvals.real)
            # plot blocked area
            sub1.fill_between((minl*maxx/10, -minl*maxx/10), miny*1.1, maxy*1.1,
                              facecolor='red', alpha=0.15)

            choice = np.logical_and(np.greater_equal(self.eigvals.real, minl),
                                    np.less_equal(self.eigvals.real, maxl))
            delayed = np.extract(choice, self.eigvals.real)
            sub2.scatter(delayed.real, delayed.imag,
                         marker='o', color='red')

            # add fundamental
            val, vect = eigenproblem.getfundamental(self, lambdas)
            sub2.scatter(0, val, marker='*', s=100, color='blue')

            # plot fundamental
            sub2.set_ylim([-1, 1])

            for la in lambdas:
                sub2.axvline(-la, color='k', linestyle='--', linewidth=0.5)

            xlo1, xup1 = sub1.get_xlim()
            ylo1, yup1 = sub1.get_ylim()
            xlo2, xup2 = sub2.get_xlim()
            ylo2, yup2 = sub2.get_ylim()
            # connection patch for first axes
            con1 = ConnectionPatch(xyA=(minl*maxx/10, (ylo1+yup1)/2), coordsA=sub1.transData,
                                   xyB=(xlo2, yup2), coordsB=sub2.transData,
                                   color='red', alpha=0.2)
            # Add left side to the figure
            fig.add_artist(con1)

            # connection patch for first axes
            con2 = ConnectionPatch(xyA=(-minl*maxx/10, (ylo1+yup1)/2), coordsA=sub1.transData,
                                   xyB=(xlo2, ylo2), coordsB=sub2.transData,
                                   color='red', alpha=0.2)
            # Add right side to the figure
            fig.add_artist(con2)

            sub2.set_xlabel('$Re(\%s)$' % label)
            sub2.set_ylabel('$Im(\%s)$' % label)
            sub2.ticklabel_format(axis='x', scilimits=[-5, 5])
            sub2.ticklabel_format(axis='y', scilimits=[-5, 5])
            if grid is True:
                sub2.grid(alpha=0.2)
            plt.tight_layout()

    def polarspectrum(self, ax=None):
        ax = ax or plt.gca()
        plt.polar(np.angle(self.eigvals[1:]), abs(self.eigvals[1:]),
                  marker='o', color='red')
        # plot fundamental
        val, vect = eigenproblem.getfundamental(self)
        plt.polar(np.angle(val), abs(val), marker='*', color='blue')

    def getfundamental(self, lambdas=None):
        if self.problem in ['kappa', 'gamma']:
            idx = 0
        elif self.problem == 'delta':
            # select real eigenvalues
            reals = self.eigvals[self.eigvals.imag == 0]
            # FIXME clean spurious big eigenvalues (temporary patch)
            reals = reals[abs(reals) < 1E3]

            if np.all(reals <= 0):
                reals = reals[reals != 0]
                fund = min(reals)
                idx = np.where(self.eigvals == fund)[0][0]
            elif np.all(reals > 0):
                fund = min(reals)
                idx = np.where(self.eigvals == fund)[0][0]
            else:
                # FIXME patch to find delta
                reals = reals[reals > 0]
                reals = reals[reals < 10]
                minreal = np.where(reals == reals.max())
                idx = np.where(self.eigvals == reals[minreal])[0][0]

        else:
            if self.problem == 'alpha':
                # select real eigenvalues
                reals = self.eigvals[self.eigvals.imag == 0]
                # all negative
                if np.all(reals < 0):
                    reals_abs = abs(self.eigvals[self.eigvals.imag == 0])
                    whichreal = np.where(reals_abs == reals_abs.min())
                else:
                    reals_abs = self.eigvals[self.eigvals.imag == 0]
                    whichreal = np.where(reals_abs == reals_abs.max())
                # blended
                idx = np.where(self.eigvals == reals[whichreal])[0][0]
            elif lambdas is not None:
                # select real eigenvalues
                choice = np.logical_and(np.greater_equal(self.eigvals.real,
                                                         min(-lambdas)),
                                        np.less_equal(self.eigvals.real,
                                                      max(-lambdas)))
                prompt = np.extract(~choice, self.eigvals)
                reals = prompt[prompt.imag == 0]
                reals_abs = abs(prompt[prompt.imag == 0])
                minreal = np.where(reals_abs == reals_abs.min())
                idx = np.where(self.eigvals == reals[minreal])[0][0]

        eigenvalue, eigenvector = self.eigvals[idx], self.eigvect[:, idx]
        return eigenvalue, eigenvector

    def get(self, geom, group, moment=0, mode=0, family=0, normalise=True,
            precursors=False):
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
        if precursors and self.problem == 'omega':
            moment = self.nA+1
            family = family-1
            iF = self.nS*family
        else:
            iF = 0

        Neven = sum(np.arange(0, moment) % 2 == 0)
        Nodd = moment-Neven

        nS = self.nS
        nE = self.nE

        iseven = 0 if moment % 2 else 1
        if Nodd >= 0:
            iS = Neven*(nE*nS)+Nodd*(nE*(nS-1))+(group-1)*((nS-1)+iseven)
        else:
            iS = group*(Neven*nE*nS)
        iS = iS+iF
        iE = iS+(nS-1)+iseven+iF

        lambdas = geom.getxs('lambda') if self.problem == 'omega' else None

        # take eigenvector
        if mode == 0:
            _, vect = eigenproblem.getfundamental(self, lambdas)
        else:
            vect = self.eigvect[:, mode]

        if normalise:
            # normalisation constant computed over total flux
            A = 1/np.linalg.norm(vect[0:nE*nS])
            vect = A*vect

        yr = np.real(vect[iS:iE])
        yi = np.imag(vect[iS:iE])
        return yr, yi
