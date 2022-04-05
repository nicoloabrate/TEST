"""
Author: N. Abrate.

File: get_energy_grid.py

Description: Utility to get default energy grid stored in
             ../material/datalib/group_structures.
"""
import numpy as np
from pathlib import Path


def get_energy_grid(grid_name):
    pwd = Path(__file__).parent.parent
    egridpath = pwd.joinpath('datalib', 'group_structures',
                             '{}.txt'.format(grid_name))
    energygrid = np.loadtxt(egridpath)
    energygrid.sort()
    return energygrid
