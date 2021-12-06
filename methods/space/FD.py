"""
Author: N. Abrate.

File: FD.py

Description: Finite Differences scheme for spatial derivates.
"""
import numpy as np


def zero(ge, f, meshtype='edges'):
    """
    Evaluate a function with Central Difference Scheme.

    Parameters
    ----------
    ge : object
        Geometry object.
    f : list[float]
        List of coefficient fues.
    meshtype : string, optional
        Mesh or staggered mesh type. The default is 'mesh'.

    Returns
    -------
    fx : np.ndarray
        1-D array with function/constant efuated over the mesh.

    """
    dicob = ge.__dict__
    N = len(dicob[meshtype])
    fx = np.zeros((1, N))
    q = 0
    NL = ge.nLayers
    if isinstance(f, (int, float, np.float)):
        f = [f]*NL
    for i in range(0, NL):
        pts = dicob[meshtype][q::]
        bord = ge.layers[i+1]
        inner_pts = np.where(pts <= bord)[0]+q*(i > 0)
        q = q+len(inner_pts)
        fx[0, inner_pts[0]:inner_pts[-1]+1] = f[i]*np.ones((1, len(inner_pts)))
        if meshtype == 'edges':
            if NL > 1 and i < NL-1 and ge.mesh[inner_pts[-1]] <= ge.layers[i+1]:
                fx[0, inner_pts[-1]] = avg(f[i], f[i+1], ge.dx[i]/2, ge.dx[i+1]/2)
    return fx


def first(ge, f, meshtype='edges', stag=True):
    """
    Evaluate first-order derivatives with Central Difference Scheme.

    Parameters
    ----------
    ge : object
        Geometry object.
    f : list[float]
        List of coefficient fues.
    meshtype : string, optional
        Mesh or staggered mesh type. The default is 'mesh'.
    stag : bool, optional
        Flag to determine where the derivative is evaluated. If ``True``, the
        derivative is evaluated on the i-th node considering the (i+1/2)-th and
        (i-1/2)-th nodes. If ``False``, the derivative is evaluated on the i-th
        node considering the (i+1)-th and (i-1)-th nodes. Default is ``True``.
    Returns
    -------
    dfdx : np.ndarray
        N-D array approximating a first-order derivative in space.

    """
    dicob = ge.__dict__
    N = len(dicob[meshtype])
    dfdx = np.zeros((2, N))
    q = 0
    NL = ge.nLayers
    mesh = dicob[meshtype]
    # ensure dimension consistency
    f = [f]*NL if isinstance(f, (int, float)) else f
    m = 1 if stag else 2

    for i in range(NL):
        pts = mesh[q::]
        dx = ge.dx
        bord = ge.layers[i+1]
        inner_pts = np.where(pts <= bord)[0]+q*(i > 0)
        P = len(inner_pts)
        q = q+P

        dfdx[0, inner_pts[0]:inner_pts[-1]+1] = -f[i]/(m*dx[i])*np.ones((P, 1)).ravel()
        # at the boundaries both right and left dxs are needed
        if stag and meshtype == 'edges':
            if NL > 1 and i < NL-1 and mesh[inner_pts[-1]] <= ge.layers[i+1]:
                dfdx[0, inner_pts[-1]] = -f[i]/(m*(dx[i]/2+dx[i+1]/2))

        dfdx[1, :] = -dfdx[0, :]

    return dfdx


def second(ge, f, meshtype='edges'):
    """
    Evaluate second-order derivatives with Central Difference Scheme.

    Parameters
    ----------
    ge : object
        Geometry object.
    f : list[float]
        List of coefficient fues.
    meshtype : string, optional
        Mesh or staggered mesh type. The default is 'mesh'.

    Returns
    -------
    d2fdx2 : np.ndarray
        N-D array approximating a second-order derivative in space.

    """
    dicob = ge.__dict__
    N = len(dicob[meshtype])
    NL = ge.nLayers
    d2fdx2 = np.zeros((3, N))
    N_old = 0
    for i in range(0, NL):
        dx = ge.dx[i]
        N = int(N_old+ge.Nx[i])
        # upper diagonal
        d2fdx2[0, N_old:N] = -f[i]/dx**2*np.ones((int(ge.Nx[i]), 1)).ravel()
        # main diagonal
        d2fdx2[1, N_old:N] = 2*f[i]/dx**2*np.ones((int(ge.Nx[i]), 1)).ravel()
        # lower diagonal
        d2fdx2[2, N_old:N] = -f[i]/dx**2*np.ones((int(ge.Nx[i]), 1)).ravel()

        if NL > 1 and i < NL-1:
            d2fdx2[2, N-1] = -f[i+1]/ge.dx[i+1]/(dx/2+ge.dx[i+1]/2)
            d2fdx2[1, N-1] = (f[i]/dx+f[i+1]/ge.dx[i+1])/(dx/2+ge.dx[i+1]/2)
            d2fdx2[0, N-1] = -f[i]/dx/(dx/2+ge.dx[i+1]/2)

        N_old = N
    # shift sub-diagonal
    d2fdx2[0, 0:-1] = d2fdx2[0, 1:]
    return d2fdx2


def avg(C1, C2, d1=1, d2=1):
    """
    Evaluate weighted averaged properties.

    Parameters
    ----------
    C1 : float
        First constant to be averaged.
    C2 : float
        Second constant to be averaged.
    d1 : float, optional
        First weight. The default is 1.
    d2 : float, optional
        Second weight. The default is 1.

    Returns
    -------
    float
        Weighted average.

    """
    return (C1*d1+C2*d2)/(d1+d2)
