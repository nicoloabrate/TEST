"""
Author: N. Abrate.

File: get_energy_grid.py

Description: Utility to get default energy grid stored in
             ../material/datalib/group_structures.
"""
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt


def get_energy_grid(grid_name):
    pwd = Path(__file__).parent.parent
    egridpath = pwd.joinpath('datalib', 'group_structures',
                             '{}.txt'.format(grid_name))
    energygrid = np.loadtxt(egridpath)
    energygrid.sort()
    return energygrid


def eplot(egrid, spectrum, ax=None, title=None, figname=None, imag=False, 
          lethargynorm=True, logx=True, logy=True, displaygrid=False,
          **kwargs, ):
    """
    Plot solution along energy for a certain portion of the phase space.

    Parameters
    ----------
    egrid : np.array or list
        Energy grid
    spectrum : np.array or list
        Spectrum plotted against energy grid
    ax : TYPE, optional
        DESCRIPTION. The default is None.
    title : TYPE, optional
        DESCRIPTION. The default is None.
    imag : TYPE, optional
        DESCRIPTION. The default is False.
    normalisation : TYPE, optional
        DESCRIPTION. The default is True.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    if len(spectrum) < len(egrid):
        spectrum = np.insert(spectrum, [0], spectrum[0])

    if logx and logy:
        loglog = True
    else:
        loglog = False

    if lethargynorm:
        u = np.log(egrid[0]/egrid)
        spectrum = spectrum/np.diff(u)

    ax = ax or plt.gca()

    plt.step(egrid, spectrum, where='pre', **kwargs)

    if loglog or logx:
        ax.set_xscale('log')
    if loglog or logy:
        ax.set_yscale('log')

    if displaygrid:
        for e in egrid:
            ax.axvline(e, c='k', lw=0.5, ls=':')

    plt.grid(which='both', alpha=0.2)
    ax.set_xlabel('E [MeV]')
    if title is not None:
        plt.title(title)

    if figname:
        plt.tight_layout()
        plt.savefig(f"{figname}.pdf")