"""
Author: N. Abrate.

File: material.py

Description: Class to handle different material regions.
"""
import numpy as np
from os import path
from pathlib import Path
from serpentTools import read
from serpentTools.settings import rc

rc['xs.reshapeScatter'] = True
rc['xs.getB1XS'] = False
rc['xs.variableGroups'] = ['kinetics', 'xs', 'xs-prod']
# This is hardcoded but depends on Serpent 2 highest Legendre polynomial expansion order
scatt_keys = [*list(map(lambda z: "infS"+str(z), range(0, 8))),
              *list(map(lambda z: "infSp"+str(z), range(0, 8)))]
xsdf_keys = ['infTot', 'infAbs', 'infDiffcoef', 'infTranspxs', 'infCapt',
             'infRemxs', 'infFiss', 'infNsf']
ene_keys = ['infNubar', 'infInvv', 'infKappa', 'infInvv',  'infChit',
            'infChip', 'infChid']

serp_keys = [*scatt_keys, *xsdf_keys, *ene_keys]

sumxs = ['Tot', 'Abs', 'Remxs']
indepdata = ['Capt', 'Fiss', 'S0', 'Nubar', 'Diffcoef', 'Chid', 'Chip']
basicdata = ['Fiss', 'Nubar', 'S0', 'Chit']
kinetics = ['lambda', 'beta', 'Chid', 'Chip']
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

            if not path.exists(path.dirname(datapath)):
                pwd = Path(__file__).parent
                fname = pwd.joinpath('datalib', '%gG' % nE, '%s' % fname)

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
            self.beta = res.resdata['fwdAnaBetaZero'][::2]
            # this to avoid confusion with python lambda function
            selfdic['lambda'] = res.resdata['fwdAnaLambda'][::2]
            self.nE = nE
            self.UniName = uniName

        except FileNotFoundError:
            try:  # read ad-hoc file format
                self._readtxt(str(fname).split('_res.m')[0], nE)
                self.nE = nE
                self.UniName = uniName
            except FileNotFoundError:
                raise OSError('Material file for %s does not exist!' % uniName)

        try:
            self.NPF = (self.beta).size
        except AttributeError:
            print('Kinetic parameters not available!')
        
        # --- complete data and perform sanity check
        L = 0
        datastr = list(self.__dict__.keys())
        # //2 since there are 'S' and 'Sp'
        S = sum('S' in s for s in datastr)//2
        self.L = S if S > L else L  # get maximum scattering order
        self.datacheck()

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

                data = np.asarray([float(val) for val in line.split()])
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
                    elif matrix.shape == (G, ):
                        selfdic[key] = matrix
                else:
                    # single-line data (scattering matrix)
                    selfdic[key] = np.asarray(data)

    def getxs(self, key, pos1=None, pos2=None):
        """
        Get material data (for a certain energy group, if needed).

        Parameters
        ----------
        key : str
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
            1-D ``numpy.ndarray`` with G/NPF (groups) rows.

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

    def perturb(self, what, howmuch, depgro=None):
        """
        
        Perturb material composition.

        Parameters
        ----------
        what : TYPE
            DESCRIPTION.
        howmuch : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        for g in range(0, self.nE):
            # no perturbation
            if howmuch[g] == 0:
               continue

            mydic = self.__dict__
            if what in indepdata:
               # update perturbed parameter
               if depgro is None:
                   delta = mydic[what][g]*howmuch[g]
                   mydic[what][g] = mydic[what][g]+delta
               else:  # select departure group for scattering matrix
                   delta = mydic[what][depgro][g]*howmuch[g]
                   mydic[what][depgro][g] = mydic[what][depgro][g]+delta

               # computesumxs = True
               # select case to ensure data consistency
               if what == 'Fiss':
                   self.Nsf[g] = self.Nubar[g]*mydic[what][g]
               elif what == 'Nubar':
                   self.Nsf[g] = self.Fiss[g]*mydic[what][g]
                   # computesumxs = False
               elif what.startswith('Chi'):
                   computesumxs = False
               elif what == 'Diffcoef':
                   # Hp: change in diffcoef implies change in capture
                   delta = 1/(3*mydic[what][g])-self.Transpxs[g]
               elif what == 'S0':
                   # change higher moments, if any
                   for l in range(0, self.L):
                       R = (mydic[what][g]/mydic[what][g]-delta)
                       key = 'S%d' % l
                       mydic[key][depgro][g] = mydic[key][depgro][g]*R
                   
               # if computesumxs is True:
               #     # force consistency
               #     for k in sumxs:
               #         if depgro != g:
               #             mydic[k][g] = mydic[k][g]+delta
               #         else:
               #             # absorption and removal indep. on in-group scatt.
               #             if k == 'Tot':
               #                 mydic[k][g] = mydic[k][g]+delta
               # # force consistency diffusion coefficient
               # sTOT = sum(self.S0[g, :])  # total in-group scattering
               # self.Transpxs[g] = self.Tot[g]-self.mu0[g]*sTOT
               # self.Diffcoef[g] = 1/(3*self.Transpxs[g])
               # # compute secondaries per collision
               # self.secpercoll[g] = (sTOT+self.Nsf[g])/(self.Tot[g])

            else:
               raise OSError('%s cannot be perturbed directly!' % what)
   
    def datacheck(self):
        """
        Check data consistency and add missing data.

        Returns
        -------
        None.

        """
        datadic = self.__dict__
        datavail = datadic.keys()
        # check basic reactions existence
        for s in basicdata:
            if s not in datavail:
                raise OSError('%s is missing in %s data!' % (s, self.UniName))
        # --- compute in-group scattering
        InScatt = np.diag(self.S0)
        sTOT = self.S0.sum(axis=1) if len(self.S0.shape) > 1 else self.S0
        # --- compute fission production cross section
        self.Nsf = self.Fiss*self.Nubar
        # --- compute missing sum reactions
        if 'Capt' in datavail:
            self.Abs = self.Fiss+self.Capt
        elif 'Abs' in datavail:
            self.Capt = self.Abs-self.Fiss
        elif 'Tot' in datavail:
            self.Capt = self.Tot-sTOT-self.Fiss    
            self.Abs = self.Fiss+self.Capt

        self.Remxs = self.Abs+sTOT-InScatt
        self.Tot = self.Remxs+InScatt
        # --- compute secondaries per collision
        self.secpercoll = (sTOT+self.Nsf)/(self.Tot)
        # --- compute mean of scattering cosine
        if 'S1' in datavail:
            self.mu0 = self.S1/self.S0
        elif 'Diffcoef' in datavail:
            tmp = 1/(3*self.Diffcoef)
            self.mu0 = (self.Tot-tmp)/sTOT
        else:
            self.mu0 = np.zeros((self.nE, ))
        # check consistency
        if abs(self.mu0).max() > 1:
            raise OSError('Average cosine larger than 1! Check input data!')
        # --- compute transport xs
        self.Transpxs = self.Tot-self.mu0*sTOT
        # --- compute diffusion coefficient
        self.Diffcoef = 1/(3*self.Transpxs)
        # --- compute diffusion length
        self.DiffLength = np.sqrt(self.Diffcoef/self.Remxs)
        # --- compute mean free path
        self.MeanFreePath = 1/self.Tot
        # --- ensure consistency kinetic parameters (if fissile medium)
        if self.Fiss.max() > 0:
            if abs(self.Chit.sum() - 1) > 1E-5:
                raise OSError('Total fission spectra in %s not normalised!' % self.UniName)
    
            for s in kinetics:
                if s in datavail:
                    kincons = True
                else:
                    kincons = False
    
            if kincons is True:
                if len(self.Chid.shape) == 1:
                    # each family has same emission spectrum
                    self.Chid = np.asarray([self.Chid]*self.NPF)
                elif self.Chid.shape != (self.NPF, self.nE):
                    raise OSError('Delayed fiss. spectrum should be (%d, %d)'
                                  % (self.NPF, self.nE))

                for g in range(0, self.nE):
                    chit = (1-self.beta.sum())*self.Chip[g]+np.dot(self.beta, self.Chid[:, g])
                    if abs(self.Chit[g]-chit) > 1E-5:
                        raise OSError('Fission spectra or delayed fractions in %s not consistent!' % self.UniName)
            