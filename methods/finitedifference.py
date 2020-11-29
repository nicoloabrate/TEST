"""
Author: N. Abrate.

File: finitedifference.py

Description: Functions for simplified 1D geometries.
"""
import numpy as np


def zero(obj, val, meshtype='mesh'):
    """
    Evaluate a function with Central Difference Scheme.

    Parameters
    ----------
    obj : object
        Geometry object.
    val : list[float]
        List of coefficient values.
    meshtype : string, optional
        Mesh or staggered mesh type. The default is 'mesh'.

    Returns
    -------
    fx : np.ndarray
        1-D array with function/constant evaluated over the mesh.

    """
    dicob = obj.__dict__
    N = len(dicob[meshtype])
    fx = np.zeros((1, N))
    q = 0

    for iLay in range(0, obj.nLayers):
        pts = dicob[meshtype][q::]
        bord = obj.layers[iLay+1]
        inner_pts = np.where(pts <= bord)[0]+q*(iLay > 0)
        q = q+len(inner_pts)
        fx[0, inner_pts[0]:inner_pts[-1]+1] = val[iLay]*np.ones((1, len(inner_pts)))

    return fx


def first(obj, val, meshtype='mesh'):
    """
    Evaluate first-order derivatives with Central Difference Scheme.

    Parameters
    ----------
    obj : object
        Geometry object.
    val : list[float]
        List of coefficient values.
    meshtype : string, optional
        Mesh or staggered mesh type. The default is 'mesh'.

    Returns
    -------
    dfdx : np.ndarray
        N-D array approximating a first-order derivative in space.

    """
    dicob = obj.__dict__
    N = len(dicob[meshtype])
    dfdx = np.zeros((2, N))
    q = 0

    mesh = dicob[meshtype]
    for iLay in range(0, obj.nLayers):
        pts = mesh[q::]
        bord = obj.layers[iLay+1]
        inner_pts = np.where(pts <= bord)[0]+q*(iLay > 0)
        q = q+len(inner_pts)

        dfdx[0, inner_pts[0]:inner_pts[-1]+1] = -val/obj.dx[iLay]*np.ones((len(inner_pts), 1)).ravel()

        if obj.nLayers > 1 and iLay < obj.nLayers and mesh[inner_pts[-1]] < obj.layers[iLay+1]:
            dfdx[0, -1] = -2*val/(obj.dx[iLay]+obj.dx[iLay+1])

        dfdx[1, :] = -dfdx[0, :]

    return dfdx


def second(obj, val, meshtype='mesh'):
    """
    Evaluate second-order derivatives with Central Difference Scheme.

    Parameters
    ----------
    obj : object
        Geometry object.
    val : list[float]
        List of coefficient values.
    meshtype : string, optional
        Mesh or staggered mesh type. The default is 'mesh'.

    Returns
    -------
    d2fdx2 : np.ndarray
        N-D array approximating a second-order derivative in space.

    """
    dicob = obj.__dict__
    N = len(dicob[meshtype])

    d2fdx2 = np.zeros((3, N))
    N_old = 0
    for iLay in range(0, obj.nLayers):
        N = N_old+obj.N[iLay]
        d2fdx2[0, 0+N_old:N] = val[iLay]*np.ones((obj.N[iLay], 1))

        if obj.nLayers > 1 and iLay < obj.nLayers:
            d2fdx2[0, N] = avg(val[iLay], val[iLay+1], obj.dx[iLay], obj.dx[iLay])

        N_old = N

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
