"""
Author: N. Abrate.

File: config.py

Description: Class to handle 1D geometry configurations.
"""
import json
from posixpath import split
import numpy as np
from TEST.utils import myjson
from TEST.geometry import Slab
from collections import OrderedDict


keys = ['split', 'layers', 'regions', 'BCs', 'energygrid', 'AngOrd', 
        'spatial_scheme']
kwkeys = ['datapath', 'verbose']
changekeys = ['Nx', 'dx', '_split', 'centers', 'edges', 'mesh',
              'ghostmesh', 'nS', 'nLayers', 'regions']
class Config:
    """Define configuration object."""

    def __init__(self, inp):

        if isinstance(inp, str) is False:
            raise OSError("Input file path is missing!")

        # parse .json file
        try:
            with open(inp) as f:
                try:
                    inp = json.load(f)
                except json.JSONDecodeError as err:
                    print(err.args[0])
                    raise OSError(err.args[0]+' in %s' % inp)
        except FileNotFoundError:
            raise OSError(f"File {inp} is missing!")
        # parse .json dict
        inp = inp['NE']
        
        args = {}
        for k in keys:
            if k in inp.keys():
                args[k] = inp[k]
            else:
                raise OSError('{} key missing in input file!'.format(k))

        kwargs = {}
        for k in kwkeys:
            if k in inp.keys():
                kwargs[k] = inp[k]

        # define initial geometry configuration
        NEassemblytypes = {}
        ireg = 0
        geometry = Slab(**args, **kwargs)
        # for steady state cases
        for reg, data in geometry.regions.items():
            if reg not in NEassemblytypes.values():
                NEassemblytypes[reg] = data

        self.time = [0]
        self.config = {0: geometry}
        if 'tEnd' in inp.keys():
            self.TimeEnd = inp['tEnd']
        else:
            self.TimeEnd = 0

        if 'nSnap' in inp.keys():
            self.nSnap = inp['nSnap']
            if isinstance(self.nSnap, (float, int)):
                dt = self.TimeEnd/self.nSnap
                self.TimeSnap = np.arange(0, self.TimeEnd+dt, dt) if dt > 0 else [0]
            elif isinstance(self.nSnap, list) and len(self.nSnap) > 1:
                self.TimeSnap = self.nSnap
            else:
                raise OSError('nProf in .json file must be list, float or int!')
        else:
            self.nSnap = 1
            self.TimeSnap = [0]
        # define configurations
        if 'config' in inp.keys():
            config = inp['config']
            for it, t in enumerate(config.keys()):
                # create geometry
                if t not in self.time:
                    self.time.append(float(t))
                else:
                    continue
                self.config[float(t)] = Slab(**args, **kwargs)
                if it > 0:
                    self.config[float(t)].regions = self.config[self.time[it]].regions
                if 'perturb' in config[t].keys():
                    self.config[float(t)].perturb(config[t]['perturb'])
                elif 'replace' in config[t].keys():
                    self.config[float(t)].replace(config[t]['replace'])

                # add new regions, if any
                for reg, data in self.config[float(t)].regions.items():
                    if reg not in NEassemblytypes.values():
                        NEassemblytypes[reg] = data

            # force same spatial cuts and split for each configuration
            lastconf = self.config[self.time[-1]]
            for nconf, conf in enumerate(self.config.values()):
                if nconf < len(self.config.values()):
                    for k in changekeys:
                        conf.__dict__[k] = lastconf.__dict__[k]
                    # update regions
                    regionmap, regionwhere = [], []
                    for x1, x2 in zip(lastconf.layers[:-1], lastconf.layers[1:]):
                        # span all old layer coordinates
                        for idr, (l, r) in conf.regionwhere.items():
                            if x1 >= l and x2 <= r:
                                regionmap.append(conf.regionmap[idr])
                                regionwhere.append((x1, x2))
                    # update regionmap and regionwhere
                    conf.regionmap = OrderedDict(zip(range(conf.nLayers), regionmap))
                    conf.regionwhere = OrderedDict(zip(range(conf.nLayers), regionwhere))

        self.NEassemblytypes = NEassemblytypes
        # check for other arguments
        if 'nSnap' in inp.keys():
            self.nSnap = inp['nSnap']

        if 'regionsplot' in inp.keys():
            self.regionslegendplot = inp['regionsplot']

        # compute nodes according to FRENETIC mesher (mesh1d.f90)
        split = self.config[0].Nx
        mesh0 = np.asarray(inp['layers'])
        nelz = split.sum()
        self.mesh = np.zeros((nelz+1,))
        idz = 0
        for iz in range(len(mesh0)-1):
            dz = (mesh0[iz+1]-mesh0[iz])/split[iz]
            for i in range(split[iz]+1):
                self.mesh[idz+i] = mesh0[iz]+i*dz
            idz += split[iz]
        self.nodes = (self.mesh[1:]+self.mesh[:-1])/2 