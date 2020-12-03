"""
Author: N. Abrate.

File: geometry.py

Description: Class for simplified 1D geometries.
"""
import numpy as np
from matplotlib.pyplot import gca
from TEST.material import Material
from collections import OrderedDict


class Slab:
    """Define slab geometry object."""

    def __init__(self, NMFP, layers, matlist, BCs, G, datapath=None):

        # assign "test" path for data library
        self.datapath = datapath
        # assign number of layers
        self.nLayers = len(layers)-1
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
        self.BC = BCs

        # assign material properties
        self.regions = {}
        minmfp = []

        for iLay in range(0, self.nLayers):

            uniName = matlist[iLay]

            if uniName in self.regions.keys():
                continue

            if datapath is not None:
                datapath = datapath.joinpath('%s' % uniName)

            self.regions[uniName] = Material(uniName, G, datapath=None)
            minmfp.append(1/np.max(self.regions[uniName].Tot))

        # assign mesh, ghost mesh and N
        Slab.mesher(self, minmfp)

        self.regionmap = OrderedDict(zip(range(0, len(matlist)), matlist))
        self.G = G
        self.geometry = 'slab'

    def mesher(self, minmfp):
        """
        Mesh the slab domain. A shadow-staggered grid is also generated.

        Returns
        -------
        None.

        """
        # preallocation
        dx = np.zeros((self.nLayers, 1))
        N = np.zeros((self.nLayers, 1))
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
        self.NT = int(sum(N))
        self.dx = dx
        self.mesh = grid
        self.stag_mesh = gridg

    def plot(self, ax=None, yVals=None, xlabel=None):
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
                vals = np.full((self.G, self.G, self.nLayers), None)
                # loop over regions
                for ireg, reg in self.regionmap.items():
                    vals[:, :, ireg] = self.regions[reg].getxs(key, pos1, pos2)

            elif key == 'beta' or key == 'lambda':

                for ireg, reg in self.regionmap.items():
                    NPF = self.regions[reg].NP

                vals = np.full((NPF, self.nLayers), None)
                # loop over regions
                for ireg, reg in self.regionmap.items():
                    vals[:, ireg] = self.regions[reg].getxs(key, pos1, pos2)

            else:
                vals = np.full((self.G, self.nLayers), None)
                # loop over regions
                for ireg, reg in self.regionmap.items():
                    vals[:, ireg] = self.regions[reg].getxs(key, pos1, pos2)

        else:

            if key.startswith('S') or key.startswith('Sp'):
                vals = np.full((self.G, self.G), None)
                # loop over regions
                vals = self.regions[reg].getxs(key, pos1, pos2)

            else:
                vals = np.full((self.G, 1), None)
                # loop over regions
                for reg, ireg in self.regionmap.items():
                    vals = self.regions[reg].getxs(key, pos1, pos2)

        return vals