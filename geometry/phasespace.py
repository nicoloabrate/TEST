"""
Author: N. Abrate.

File: PhaseSpace.py

Description: Class to handle phase space operations.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import ConnectionPatch
from scipy.special import eval_legendre


class PhaseSpace:
    """Define phase space object."""

    def __init__(self, geometry, solution, energygrid=None, source=False,
                 normalize=False, whichnorm='phasespace', operators=None):
        """
        Initialise phase space object.

        Parameters
        ----------
        geometry : object
            Geometry object containing the material data and the spatial mesh.
        solution : dict or ndarray
            Solution living in the phase space. If it is a dict, it is treated 
            as the solution of an eigenvalue problem. If it is an ndarray and
            ``source`` is True, it is treated as the solution of a source-driven
            problem.
        energygrid : ndarray, optional
            Multi-group energy grid structure. The default is None.
        source : bool, optional
            If ``True``, the solution is not normalised. 
            The default is ``False``.

        Returns
        -------
        None.

        """
        self.geometry = geometry

        if energygrid is not None:
            self.energygrid = energygrid

        if isinstance(solution, dict):
            if operators is not None:
                self.model = operators.model
            else:
                raise OSError('Numerical operators object is needed!')
            self.eigvals = solution['eigenvalues']
            self.eigvect = solution['eigenvectors']
            self.problem = solution['problem']

            self.nA = operators.nA
            self.nE = operators.nE
            self.nS = operators.nS

            if normalize is True:
                self.normalize(which=whichnorm)
        elif isinstance(solution, (np.ndarray)) and source is True:
            self.flux = solution
        else:
            raise OSError('Type {} cannot be handled by phase space!'.format(type(solution)))

    def braket(self, v1, v2=None):
        """
        Compute bra-ket product over the phase space.

        Parameters
        ----------
        v1 : ndarray
            DESCRIPTION.
        v2 : ndarray
            DESCRIPTION.
        geom : object
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if self.geometry.AngOrd > 0:
            raise OSError('GPT cannot be applied yet to transport solutions!')

        v1v2 = np.multiply(v1, v2) if v2 is not None else v1
            
        G = self.geometry.nE
        S = self.geometry.nS
        grid = self.geometry.mesh
        I = 0
        # TODO consider integration over non-group energy grid
        for g in range(0, G):
            skip = g*S
            I = I+np.trapz(v1v2[skip:skip+S], x=grid)
        return I

    def interp(self, yp, xx, isref=True):
        """
        Interpolate linearly vector yp on a different phase space geometrical 
        grid.

        Parameters
        ----------
        yp : ndarray
            Vector to be interpolated.
        xx : TYPE
            Old or new grid according to isref argument.
        isref : bool, optional
            Flag to choose new grid. If ``True`` the new grid is assumed
            to be the one in self.geometry. The default is ``True``.

        Raises
        ------
        OSError
            Phase space interpolation cannot be applied yet to transport solutions!

        Returns
        -------
        y : ndarray
            Vector interpolated over the new phase space grid.

        """
        G = self.geometry.nE
        if self.geometry.AngOrd > 0:
            raise OSError('Phase space interpolation cannot be applied yet to transport solutions!')

        if isref is True:
            x = self.geometry.mesh
            xp = xx            
        else:
            x = xx
            xp = self.geometry.mesh

        n = len(xp)
        N = len(x)
        y = np.zeros((G*N,), dtype=complex)
        for g in range(0, G):
            y[g*N:(g+1)*N] = np.interp(x, xp, yp[g*n:(g+1)*n])

        return y

    def normalize(self, which='phasespace', adjoint=None, power=None,
                  A=None):
        """
        Normalize eigenvectors according to a user-defined criterion.

        Parameters
        ----------
        which : str, optional
            Normalisation condition. The default is ``None``.
        adjoint : object, optional
            Adjoint eigenfunctions for bi-orthogonal normalisation.
            The default is ``None``.

        Returns
        -------
        None.

        """
        # TODO: normalise according to numerical method employed (PN, SN, diff)
        # use getter methods to normalise over total flux
        if which == 'phasespace':
            # normalise on the inner product over the phase space
            for iv, v in enumerate(self.eigvect.T):
                C = self.braket(self.eigvect[:, iv], self.eigvect[:, iv])
                self.eigvect[:, iv] = v/np.sqrt(C)
        elif which == 'norm2':
            # normalise to have unitary euclidean norm
            for iv, v in enumerate(self.eigvect.T):
                self.eigvect[:, iv] = v/np.linalg.norm(v)
        elif which == 'power':
            # normalise to have fixed power
            power = 1 if power is None else power
            print('Under development')
        # elif which == 'totalflux':
        #     # normalise total flux
        #     if self.model == 'Diffusion':

        #     elif self.model == 'PN':
        #         self.get()
        #     elif self.model == 'SN':

        # elif which == 'peaktotalflux':

        elif which == 'reaction':
            # normalise to have unitary reaction rate
            for iv, v in enumerate(self.eigvect.T):
                self.eigvect[:, iv] = v/self.braket(A.dot(v))
        elif which == 'biorthogonal':
            # normalise to have <phix, F*phi>=<FX*phix, phi>=1
            n, m = self.solution.eigvect.shape[1], self.eigvect.shape[1]
            # check consistency
            if adjoint is None:
                raise OSError('For biorthogonal normalisation the adjoint input argument is needed!')
            elif n != m:
                raise OSError('Adjoint modes number does not match forward modes! {}!={}'.format(n, m))
            elif adjoint.operators.state != 'steady':
                raise OSError('Bi-orthogonal normalisation only available for steady state condition!')

            for iv, v in enumerate(self.eigvect.T):
                v_adj = self.solution.eigvect[:, iv]
                self.eigvect[:, iv] = v/self.braket(self.operators.F.dot(v_adj), v)

    def plot(self, group, angle=0, mode=0, family=0, precursors=False,
             ax=None, title=None, imag=False, **kwargs):

        yr, yi = self.get(group, angle=angle, mode=mode, family=family,
                          precursors=precursors)
        x = self.geometry.mesh if len(yr) == self.geometry.nS else self.geometry.stag_mesh
        ax = ax or plt.gca()

        y = yr if imag is False else yi
        plt.plot(x, y, **kwargs)

        ax.locator_params(nbins=8)
        ax.set_ylim(y.min())
        ax.set_xlim([min(self.geometry.layers), max(self.geometry.layers)])
        ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

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

        val, vect = self.getfundamental(lambdas)
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
            val, vect = self.getfundamental(lambdas)
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
        val, vect = self.getfundamental()
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

    def get(self, group=None, angle=None, moment=0, mode=0, family=0,
            normalise=True, precursors=False):
        """
        Get spatial flux distribution for group, angle/moment and spatial mode.

        Parameters
        ----------
        group : int, optional
            Group number. If ``None``, all groups are displayed. Default is 
            ``None``.
        angle : int or float, optional
            Angle to be displayed. If ``int``, the direction matching the 
            number will be provided, if ``float`` the direction closest to 
            this number will be provided. Default is None.
        moment : int
            Flux moment number. Default is 0.
        mode : int
            Eigen-mode number. Default is 0.
        family : int, optional
            Precursors family number. Default is 0.
        precursors : bool, optional
            Flag for precursor concentration. Default is ``False``.

        Returns
        -------
        yr : ndarray
            Real part of the flux mode.
        yi : ndarray
            Imaginary part of the flux mode.

        """
        if self.model == 'PN' or self.model == 'Diffusion':
            yr, yi = self._getPN(group, angle=angle, moment=moment, mode=mode,
                                 family=family, normalise=normalise,
                                 precursors=precursors)
        elif self.model == 'SN':
            yr, yi = self._getSN(group, angle=angle, moment=moment, mode=mode,
                                 family=family, normalise=normalise,
                                 precursors=precursors)
        return yr, yi

    def _getPN(self, group, angle=None, moment=0, mode=0, family=0,
               normalise=True, precursors=False):
        """
        Get spatial flux distribution for group, angle and spatial mode.

        Parameters
        ----------
        geom : object
            Geometry object.
        angle : int
            Moment/direction number.
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
        nE, nA, nS = self.nE, self.nA, self.nS

        if angle is None:
            if precursors and self.problem == 'omega':
                moment = nA
                nF = self.nF
                family = family-1
                iF = nS*family
            else:
                iF = 0

            No = (nA+1)//2 if nA % 2 != 0 else nA//2
            Ne = nA+1-No
    
            lambdas = self.geometry.getxs('lambda') if self.problem == 'omega' else None
    
            # take eigenvector
            if mode == 0:
                _, vect = self.getfundamental(lambdas)
            else:
                vect = self.eigvect[:, mode]
    
            # if normalise:
            #     # normalisation constant computed over total flux
            #     totflx = np.zeros((nE*nS,), dtype=complex)
            #     for gro in range(0, nE):
            #         skip = (Ne*nS+No*(nS-1))*gro
            #         totflx[gro*nS:(gro+1)*nS] = vect[skip:skip+nS]
            #     vect = vect/np.linalg.norm(totflx)
            G = self.geometry.nE
            dim = self.geometry.nS if moment % 2 == 0 else self.geometry.nS-1
            # preallocation for group-wise moments
            yr, yi = np.zeros((dim*G, )), np.zeros((dim*G, ))
            gro = np.arange(0, self.geometry.nE) if group is None else group

            for g in gro:

                if precursors is False:
                    # compute No and Ne for the requested moment/angle
                    skip = (Ne*nS+No*(nS-1))*(g-1)
                    NO = (moment+1)//2 if (moment-1) % 2 != 0 else moment//2
                    NE = moment-NO
                    M = nS if moment % 2 == 0 else nS-1
                    iS = skip+NE*nS+NO*(nS-1)
                    iE = skip+NE*nS+NO*(nS-1)+M
    
                else:
                    iS = (Ne*nS+No*(nS-1))*nE+g*nF+iF
                    iE = (Ne*nS+No*(nS-1))*nE+g*nF+iF+nS
                # store slices
                yr[g*dim:(g+1)*dim] = np.real(vect[iS:iE])
                yi[g*dim:(g+1)*dim] = np.imag(vect[iS:iE])

        else:
            # build angular flux and evaluate in angle
            tmp = np.zeros((dim*G))
            for n in range(0, self.nA):
                # interpolate to have consistent PN moments
                y = vect if n % 2 == 0 else self.interp(vect, self.geometry.stag_mesh)
                tmp = tmp+(2*n+1)/2*eval_legendre(n, angle)*y
            # get values for requested groups
            iS, iE = 0, -1 if group is None else nS*(group-1), nS*(group-1)+nS
            yr, yi = np.real(y[iS:iE]), np.imag(y[iS:iE])

        return yr, yi

    def _getSN(self, group, angle=None, moment=0, mode=0, family=0,
               normalise=True, precursors=False):
        """
        Get spatial flux distribution for group, angle and spatial mode.

        Parameters
        ----------
        geom : object
            Geometry object.
        angle : int
            Moment/direction number.
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
        nE, nA = self.nE, self.nA
        nS = self.nS if self.geometry.spatial_scheme == 'FD' else self.nS-1
        mu = self.geometry.QW['mu']
        w = self.geometry.QW['w']
        if angle is not None:

            if isinstance(angle, int):
                idx = angle
            elif isinstance(angle, float):
                # look for closest direction
                idx = np.argmin(abs(mu-angle))

            if precursors and self.problem == 'omega':
                moment = nA
                nF = self.nF
                family = family-1
                iF = nS*family
            else:
                iF = 0

            lambdas = self.geometry.getxs('lambda') if self.problem == 'omega' else None
    
            # take eigenvector
            if mode == 0:
                _, vect = self.getfundamental(lambdas)
            else:
                vect = self.eigvect[:, mode]

            G = self.geometry.nE
            # preallocation for group-wise moments
            yr, yi = np.zeros((nS*G, )), np.zeros((nS*G, ))
            gro = np.arange(0, self.geometry.nE) if group is None else group

            for g in gro:

                if precursors is False:
                    # compute No and Ne for the requested moment/angle
                    iS = (nS*idx)*(g-1)
                    iE = (nS*idx)*(g-1)+nS
                else:
                    iS = len(mu)*nS*nE+g*nF+iF
                    iE = len(mu)*nS*nE+g*nF+iF+nS
                # store slices
                yr[g*nS:(g+1)*nS] = np.real(vect[iS:iE])
                yi[g*nS:(g+1)*nS] = np.imag(vect[iS:iE])

        else:
            # build angular flux and evaluate in angle
            tmp = np.zeros((dim*G))


        return yr, yi