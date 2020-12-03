"""
Author: N. Abrate.

File: material.py

Description: Class to handle different material regions.
"""
import numpy as np
from pathlib import Path
from serpentTools import read
from serpentTools.settings import rc

rc['xs.reshapeScatter'] = True
rc['xs.getB1XS'] = False
rc['xs.variableGroups'] = ['kinetics', 'xs', 'xs-prod']

scatt_keys = [*list(map(lambda z: "infS"+str(z), range(0, 8))),
              *list(map(lambda z: "infSp"+str(z), range(0, 8)))]
xsdf_keys = ['infTot', 'infAbs', 'infDiffcoef', 'infTranspxs', 'infCapt',
             'infRemxs', 'infFiss', 'infNsf']
ene_keys = ['infNubar', 'infInvv', 'infKappa', 'infInvv',  'infChit',
            'infChip', 'infChid']

serp_keys = [*scatt_keys, *xsdf_keys, *ene_keys]


class Material():
    """Create material regions with multi-group constants."""

    def __init__(self, uniName, nE, datapath=None):
        """


        Parameters
        ----------
        uniName : TYPE
            DESCRIPTION.
        nE : TYPE
            DESCRIPTION.
        datapath : TYPE, optional
            DESCRIPTION. The default is None.

        Raises
        ------
        OSError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if datapath is None:

            pwd = Path(__file__).parent
            datapath = pwd.joinpath('datalib', '%gG' % nE,
                                    '%s_res.m' % uniName)

        try:  # read Serpent 2 output with serpentTools

            if str(datapath).endswith('_res.m'):
                fname = str(datapath)

            else:
                fname = '%s_res.m' % str(datapath)

            res = read(str(fname))
            univ = []
            for key in sorted(res.universes):
                univ.append(key[0])

            if uniName not in univ:
                raise OSError('%s does not contain %s GC_UNIVERSE!' % (datapath, uniName))

            else:
                data = res.getUniv(uniName, 0)

            if len(data.infExp['infAbs']) != nE:
                raise OSError('%s no. energy groups not match with input G!' % datapath)

            selfdic = self.__dict__

            for my_key in serp_keys:

                if my_key.startswith('infS') or my_key.startswith('infSp'):
                    vals = np.reshape(data.infExp[my_key], (nE, nE))

                else:
                    vals = data.infExp[my_key]

                selfdic[my_key.split('inf')[1]] = vals

            # kinetics parameters
            selfdic['beta'] = res.resdata['fwdAnaBetaZero'][::2]
            selfdic['lambda'] = res.resdata['fwdAnaLambda'][::2]

        except FileNotFoundError:

            try:  # read ad-hoc file format
                self._readtxt(str(fname).split('_res.m')[0], nE)

            except FileNotFoundError:
                raise OSError('Material file for %s does not exist!' % uniName)

        self.NP = (self.beta).size

    def _readtxt(self, fname, nE):
        """
        Parse the material data from a .txt file.

        Macro-group constants are parsed from a formatted file with column-wise data
        separated by headers beginning with "#" and the name of the data:
            * Tot: total cross section [cm^-1]
            * Transpxs: transport cross section [cm^-1]
                        It is defined as total_xs-avg_direction*scattering_xs
                        according to P1 approximation.
            * Diffcoef: diffusion coefficient [cm]
                        It is defined as 1/(3*Transpxs).
            * Abs: absorption cross section [cm^-1]
                   It is the sum of Capt and Fiss cross sections.
            * Capt: capture cross section [cm^-1]
            * Fiss: fission cross section [cm^-1]
            * Remxs: removal cross section [cm^-1]
                    It is the sum of Abs and group-removal.
            * Chit: total emission spectrum [-]
            * Chip: prompt emission spectrum [-]
            * Chid: delayed emission spectrum [-]
            * Nsf: fission production cross section [cm^-1]
            * Nubar: neutron multiplicities [-]
            * Kappa: average fission deposited heat [MeV]
            * Invv: particle inverse velocity [s/cm]
            * S0, S1, S2,... : scattering matrix cross section [cm^-1]
            * Sp0, Sp1, Sp2,... : scattering production matrix cross section [cm^-1]
            * beta: delayed neutron fractions [-]
            * lambda: precursors families decay constant [-]

        Parameters
        ----------
        fname : string
            Material data file name.
        nE : int
            Number of energy groups.

        Returns
        -------
        None.

        """
        selfdic = self.__dict__
        G = None

        lines = open(fname).read().split('\n')

        for il, line in enumerate(lines):

            if line.startswith('#'):
                key = (line.split('#')[1]).strip()
                matrix = None

            elif line == '':
                continue

            else:

                data = np.asarray([float(val) for val in line.split(' ')])
                if G is None:
                    G = len(data)

                if G != nE:
                    raise OSError('Number of groups in line %g is not consistent!', il)

                if key.startswith('S') or key.startswith('Sp'):
                    # multi-line data (scattering matrix)
                    if matrix is None:
                        matrix = np.asarray(data)
                    else:
                        matrix = np.c_[matrix, data]

                    if matrix.shape == (G, G):
                        selfdic[key] = matrix.T

                else:
                    # single-line data (scattering matrix)
                    selfdic[key] = np.asarray(data)

        # construct missing data, if any
        if 'Nsf' not in selfdic.keys():
            selfdic['Nsf'] = selfdic['Nubar']*selfdic['Fiss']
        if 'Abs' not in selfdic.keys():
            selfdic['Abs'] = selfdic['Fiss']+selfdic['Capt']
        if 'Tot' not in selfdic.keys():
            # evaluate total scattering XS
            totscatt = np.zeros((selfdic['Fiss'].shape))
            order = 0
            for keys in selfdic.keys():
                if key.startswith('S'):
                    try:
                        totscatt = totscatt+np.sum(selfdic['S%d' % order],
                                                   axis=0)
                        order = order+1
                    except KeyError:
                        pass

            selfdic['Tot'] = selfdic['Abs']+totscatt

    def getxs(self, key, pos1=None, pos2=None):
        """
        Get material data for a certain energy group.

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

        Returns
        -------
        vals : numpy.ndarray
            1-D ``numpy.ndarray`` with G (groups) rows.

        """
        if pos1 is None and pos2 is None:
            vals = self.__dict__[key]

        else:
            if key.startswith('S') or key.startswith('Sp'):

                if pos2 is None:
                    raise OSError('Two coordinates needed for %s data' % key)

                else:
                    vals = self.__dict__[key][pos1, pos2]

            else:
                vals = self.__dict__[key][pos1]

        return vals
