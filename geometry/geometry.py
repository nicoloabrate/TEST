"""
Author: N. Abrate.

File: geometry.py

Description: Class for simplified 1D geometries.
"""
import json
import numpy as np
from matplotlib.pyplot import gca
from TEST.material import Material
from collections import OrderedDict
from copy import deepcopy as cp
from scipy.special import roots_legendre, eval_legendre


class Slab:
    """Define slab geometry object."""

    def __init__(self, NMFP, layers, matlist, BCs, G, AngOrd, spatial_scheme,
                 datapath=None):

        # assign number of layers
        self.nLayers = len(layers)-1
        # check input consistency
        if self.nLayers != len(matlist):
            raise OSError('{} regions but only {} materials specified!'.format(self.nLayers, len(matlist)))

        # assign "test" path for data library
        self.datapath = datapath

        # assign layers coordinates
        self.layers = layers

        # set number of mean free path
        if isinstance(NMFP, (int, float)):
            self._NMFP = NMFP*np.ones((self.nLayers))

        elif isinstance(NMFP, (list, np.ndarray)):

            if len(NMFP) == self.nLayers:
                self._NMFP = NMFP

            else:
                raise OSError('NMFP specified are not consistent with' +
                              ' the number of regions. Please provide one'
                              + ' NMFP or as many NMFP as the number of' +
                              'layers.')
        else:
            raise TypeError('Cannot read type %s for NMFP!' % type(NMFP))

        # set Boundary Conditions
        self.BC = BCs if isinstance(BCs, list) else [BCs]

        # assign material properties
        self.regions = {}
        minmfp = np.zeros((self.nLayers, ))

        for iLay in range(0, self.nLayers):

            uniName = matlist[iLay]

            if uniName not in self.regions.keys():

                if datapath is not None:
                    path = datapath[iLay]
                else:
                    path = None
    
                self.regions[uniName] = Material(uniName, G, datapath=path)

            # consistency check precursor families
            if 'NPF' in self.__dict__.keys():
                if iLay == 0:
                    self.NPF = self.regions[uniName].NPF
                else:
                    if self.NPF != self.regions[uniName].NPF:
                        raise OSError('Number of precursor families in %s not' +
                                      ' consistent with other regions' % uniName)
            if AngOrd > 0:
                minmfp[iLay] = min(self.regions[uniName].MeanFreePath)
            else:
                minmfp[iLay] = np.min(self.regions[uniName].DiffLength)
        # assign mesh, ghost mesh and N
        self.mesher(minmfp)

        self.regionmap = OrderedDict(zip(range(0, len(matlist)), matlist))
        self.AngOrd = AngOrd
        self.spatial_scheme = spatial_scheme
        self.nE = G
        self.geometry = 'slab'

    def mesher(self, minmfp):
        """
        Mesh the slab domain. A shadow-staggered grid is also generated.

        Returns
        -------
        None.

        """
        # preallocation
        dx = np.zeros((self.nLayers,))
        N = np.zeros((self.nLayers,))
        old_grid = np.empty(0)
        old_grid_stag = np.empty(0)

        for iLay in range(0, self.nLayers):
            # compute grid spacing
            if max(self._NMFP) < 0:  # assign user-defined number of points
                deltalay = self.layers[iLay+1]-self.layers[iLay]
                dx[iLay] = deltalay/abs(self._NMFP[iLay])
                N[iLay] = int(np.ceil(deltalay/dx[iLay]))

            elif (np.round(self._NMFP) != self._NMFP).any():  # user-defined dx
                deltalay = self.layers[iLay+1]-self.layers[iLay]
                dx[iLay] = self._NMFP
                N[iLay] = int(np.ceil(deltalay/dx[iLay]))

            else:
                deltalay = self.layers[iLay+1]-self.layers[iLay]
                dx[iLay] = minmfp[iLay]/self._NMFP[iLay]
                N[iLay] = int(np.ceil(deltalay/dx[iLay]))

            # grid
            ngrid = np.linspace(self.layers[iLay], self.layers[iLay+1],
                                int(N[iLay]))
            grid = np.concatenate((old_grid, ngrid))
            old_grid = grid
            # ghost grid
            dx[iLay] = ngrid[1]-ngrid[0]
            ngridg = np.arange(self.layers[iLay]+dx[iLay]/2,
                               self.layers[iLay+1], dx[iLay])
            gridg = np.concatenate((old_grid_stag, ngridg))
            old_grid_stag = gridg

            if len(np.unique(grid)) < len(grid):
                N[iLay] = N[iLay]-1
                grid = np.unique(grid)

        if np.any(grid) or np.any(gridg):
            grid[grid == 0] = np.finfo(float).eps
            gridg[gridg == 0] = np.finfo(float).eps

        self.N = N
        self.nS = int(sum(N))
        self.dx = dx
        self.mesh = grid
        self.stag_mesh = gridg

    def plotmesh(self, ax=None, yVals=None, xlabel=None):
        """Plot mesh grid."""
        if yVals is None:
            ymin, ymax = 0, 1

        else:
            ymin, ymax = np.min(yVals), np.max(yVals)

        ax = ax or gca()
        ax.vlines(self.mesh, ymin, ymax, colors='k',
                  linestyles='solid')

        ax.vlines(self.stag_mesh, ymin, ymax, colors='k',
                  linestyles='dashed')

        xlabel = xlabel if xlabel is not None else 'x coordinate [cm]'
        ax.set_xlabel(xlabel)
        ax.set_xticks(self.layers)

    def displaygeom(self, ax=None, xlabel=None, labels=None):
        """Plot regions."""

        ax = ax or gca()
        c = ['royalblue', 'firebrick',
             'forestgreen', 'gold', 'darkorange',
             'darkviolet']

        c = dict(zip(self.regions.keys(), c))
        for i in range(0, self.nLayers):
            which = self.regionmap[i]
            col = c[which]
            ax.axvspan(self.layers[i], self.layers[i+1],
                       alpha=0.5, color=col, label=which)

        xlabel = xlabel if xlabel is not None else 'x coordinate [cm]'
        ax.set_xlabel(xlabel)
        ax.set_xticks(self.layers)
        if labels is None:
            ax.legend(bbox_to_anchor=(1.05, 1))
        else:
            ax.legend(labels, bbox_to_anchor=(1.05, 1))

    def getxs(self, key, pos1=None, pos2=None, region=None):
        """
        Get material data for a certain region and energy group.

        Parameters
        ----------
        key : string
            User selected nuclear data.
        pos1 : int, optional
            Departure energy group for scattering matrix. If not provided,
            data over all the energy groups are returned.
            The default is ``None``.
        pos2 : int, optional
            Arrival energy group for scattering matrix. If not provided,
            data over all the energy groups are returned.
            The default is ``None``.
        region : string, optional
            Region name. The default is ``None``.

        Returns
        -------
        vals : numpy.ndarray
            ``numpy.ndarray`` with G (groups) rows and R (regions) columns.

        """
        if region is None:

            if key.startswith('S') or key.startswith('Sp'):
                if key == 'S' or key == 'Sp':  # take all moments
                    # look for maximum number of moments available in the data
                    L = 0
                    for ireg, reg in self.regionmap.items():
                        S = self.regions[reg].L
                        L = S if S > L else L  # get maximum scattering order
                    L = 1 if L == 0 else L
                    shape = (self.nE, self.nE, self.nLayers, L)
                    allL = True
                else:
                    shape = (self.nE, self.nE, self.nLayers)
                    allL = False

                vals = np.full(shape, None)
                # loop over regions
                for ireg, reg in self.regionmap.items():
                    if allL is True:
                        old_key = 'S'
                        # get all scattering order matrices for each region
                        for l in range(0, L):
                            key = '%s%d' % (old_key, l)
                            vals[:, :, ireg, l] = self.regions[reg].getxs(key, pos1, pos2)
                    else:
                        vals[:, :, ireg] = self.regions[reg].getxs(key, pos1, pos2)

            elif key == 'beta' or key == 'lambda':

                for ireg, reg in self.regionmap.items():
                    NPF = self.regions[reg].NPF

                vals = np.full((NPF, self.nLayers), None)
                # loop over regions
                for ireg, reg in self.regionmap.items():
                    vals[:, ireg] = self.regions[reg].getxs(key, pos1, pos2)

            else:
                vals = np.full((self.nE, self.nLayers), None)
                # loop over regions
                for ireg, reg in self.regionmap.items():
                    vals[:, ireg] = self.regions[reg].getxs(key, pos1, pos2)

        else:

            if key.startswith('S') or key.startswith('Sp'):
                vals = np.full((self.nE, self.nE), None)
                # loop over regions
                vals = self.regions[reg].getxs(key, pos1, pos2)

            else:
                vals = np.full((self.nE, 1), None)
                # loop over regions
                for reg, ireg in self.regionmap.items():
                    vals = self.regions[reg].getxs(key, pos1, pos2)

        return vals

    def perturb(self, perturbation):
        """
        Add perturbations to an unperturbed, reference system.

        Parameters
        ----------
        what : str
            DESCRIPTION.
        where : ndarray
            DESCRIPTION.
        howmuch : float
            DESCRIPTION.
        system : object
            DESCRIPTION.

        Returns
        -------
        None.

        """
        # parse file content
        if isinstance(perturbation, str):
            if '.json' not in perturbation:
                perturbation = '%s.json' % perturbation
            with open(perturbation) as f:
                perturbation = json.load(f)

        iP = 0  # perturbation counter
        regs = list(self.regionmap.values())
        for k, v in perturbation.items():
            # check perturbation consistency
            if 'where' not in v.keys():
                raise OSError('Perturbation coordinates are missing')
            if 'howmuch' not in v.keys():
                raise OSError('Perturbation intensities are missing')
            if 'depgro' not in v.keys():
                v['depgro'] = None

            x = v['where']
            hw = v['howmuch']
            dg = v['depgro']
            if len(hw) != self.nE:
                raise OSError('The perturbation intensities required should be %d' % self.nE)

            if isinstance(x, tuple):
                x = [x]

            for x1, x2 in x:  # loop over perturbation coordinates
                coo = zip(self.layers[:-1], self.layers[1:])
                iP = iP + 1
                for r, l in coo:  # loop to identify regions involved
                    # perturbation inside region
                    if x1 >= r and x2 <= l:
                        # is included in this region
                        if x1 not in self.layers:
                            idx = self.layers.index(r)
                            self.layers.insert(idx+1, x1)
                        if x2 not in self.layers:
                            self.layers.insert(idx+2, x2)
                        # add new region
                        if x1 == r and x2 < l:  # on the right
                            oldreg = regs[idx]
                            regs.insert(idx, 'Perturbation%d' % iP)
                        elif (x1, x2) == (r, l):
                            idx = self.layers.index(r)
                            oldreg = regs[idx]
                            regs[idx] =  'Perturbation%d' % iP
                        elif x2 == l:  # on the left
                            oldreg = regs[idx]
                            regs.insert(idx+1, 'Perturbation%d' % iP)
                        else:
                            regs.insert(idx, regs[idx])
                            oldreg = regs[idx]
                            regs.insert(idx+1, 'Perturbation%d' % iP)

                        self.regions['Perturbation%d' % iP] = cp(self.regions[oldreg])
                        self.regions['Perturbation%d' % iP].perturb(k, hw, dg)
                    # perturbation between two or more regions
                    elif x1 < l and x2 > l:
                        raise OSError('Perturbations can be applied one region at a time!')

            self.nLayers = len(self.layers)-1
            self.regionmap = OrderedDict(zip(range(0, self.nLayers), regs))
            # update mesh to take into account new layers
            if list(set(self._NMFP.tolist())) == self._NMFP.tolist():
                self._NMFP = self._NMFP*np.ones((self.nLayers))
            else:
                raise OSError('Number of MFP must be the same for all regions for perturbation!')

            minmfp = np.zeros((self.nLayers, ))
    
            for iLay in range(0, self.nLayers):
                uniName = self.regionmap[iLay]
                if self.AngOrd > 0:
                    minmfp[iLay] = min(self.regions[uniName].MeanFreePath)
                else:
                    minmfp[iLay] = np.min(self.regions[uniName].DiffLength)

            self.mesher(minmfp)

    def computeQW(self):
        """
        Compute the set of quadrature weights and normalisation coefficients
        needed for the PN approximation.

        Returns
        -------
        None.

        """
        # compute Legendre expansion coefficients for SN
        sm = self.getxs('%s' % 'S')
        L = sm.shape[3]
        mu, w = roots_legendre(self.AngOrd)
        # ensure positive and then negative directions
        mu[::-1].sort()

        PL = np.zeros((L, self.AngOrd))
        for order in range(0, L):
            PL[order, :] = eval_legendre(order, mu)
        C = (2*np.arange(0, self.AngOrd)+1)/2
        QW = {'L': L, 'mu': mu, 'w': w, 'PL': PL, 'C': C}
        self.QW = QW