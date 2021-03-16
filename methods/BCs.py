"""
Created on Sun Mar 14 17:04:03 2021

author: N. Abrate.

file: .py

description:
"""
import os
import warnings
import numpy as np
from scipy.special import eval_legendre, roots_legendre
from sympy import Symbol, legendre, integrate

warnings.simplefilter('ignore')


def imposeBC(op, slab):
    """
    Impose boundary conditions specified by the user.

    Parameters
    ----------
    op : object
        Transport equation operators in PN approximation.
    slab : object
        Geometry object representing a 1D cartesian domain.

    Raises
    ------
    OSError
        Unknown boundary condition, if the option is not recognised

    Returns
    -------
    None.

    """
    # FIXME: to correctly handle one-side BC, only positive/negative directions
    # have to be kept!
    M = op.nS
    N = op.nA
    G = op.nE
    BCs = slab.BC
    op.BC = BCs

    # copy leakage operator to new variable
    L = op.Linf  # .tocsr()+0

    for bc in BCs:

        # FIXME: actually only the same bc can be handled on two boundaries
        # TODO: check boundary conditions consistency (if different can be imposed)


        if bc in ['zero', 'zeroflux']:
            print('Under development')
            L[0, 0] = 1
            L[0, slab.NT] = 0
            L[slab.NT-1, slab.NT-1] = 1
            L[slab.NT-1, -1] = 0

            op.F[0, 0] = 0
            op.F[slab.NT, slab.NT] = 0
            op.R[0, 0] = 0
            op.R[slab.NT, slab.NT] = 0
            op.S[0, 0] = 0
            op.S[slab.NT, slab.NT] = 0

        else:

            if bc in ['markeven', 'Markeven', 'markEven', 'MarkEven']:
                A = MarkCoeffs(N, even=True)

            elif bc in ['mark', 'Mark']:
                A = MarkCoeffs(N)

            elif bc in ['marshak', 'Marshak']:
                A = MarshakCoeffs(N)

            else:
                raise OSError('Unknown boundary condition %s!' % bc)

            A = _getcoeffs(A)
            m, n = A.shape

            A[0:m//2, :] = -2/slab.dx[0]*A[0:m//2, :]
            A[m//2:m, :] = 2/slab.dx[-1]*A[m//2:m, :]

            count = 0

            for gro in range(0, op.nE):

                for moment in range(0, n):

                    Neq = 2*moment  # only even moments need BCs
                    No = (N+1)//2 if N % 2 != 0 else N//2
                    Ne = N+1-No
                    ip = moment*(M+(M-1))+gro*(2*M-1)
                    igr = gro*(2*M-1)
                    igc = gro*(M-1)
                    ig = gro*(No*(M-1)+Ne*M)

                    if moment >= 1:  # *2 for one-side f.d. FIXME >=
                        # right boundary, lower diag
                        L[igr+ip, igc+ip-(M-1)] = 2*L[igr+ip, igc+ip-(M-1)]
                        # left boundary, lower diag
                        L[M+igr+ip-1, igc+ip-1] = 2*L[M+igr+ip-1, igc+ip-1]

                    if moment < n-1 or N % 2 != 0:  # no last and odd eq.
                        # right boundary, upper diag
                        L[igr+ip, igc+ip+(gro+1)*M] = 2*L[igr+ip, igc+ip+(gro+1)*M]
                        # left boundary, upper diag
                        L[M+igr+ip-1, (M-1)+igc+ip-1+(gro+1)*M] = 2*L[M+igr+ip-1, (M-1)+igc+ip-1+(gro+1)*M]

                    # set non-diagonal entries (even moments)
                    jj = np.arange(0, n)
                    iEv = jj*(2*M-1)

                    if moment == 0:  # 1st row, eqs 1 and 2 (Upper)
                        # right boundary, lower diag
                        L[ig, ig+iEv] = (Neq+1)/(2*Neq+1)*A[0, jj]  # angle>0
                        # left boundary, lower diag
                        L[M+ig-1, ig+iEv+M-1] = (Neq+1)/(2*Neq+1)*A[m//2, jj]  # angle<0

                    else:
                        # sum coeffs in previous row (Lower)
                        L[ip+igr, ig+iEv] = Neq/(2*Neq+1)*A[count, jj]  # angle>0
                        L[ip+igr+M-1, ig+iEv+M-1] = Neq/(2*Neq+1)*A[count+m//2, jj]  # angle<0

                        if moment < n-1 or N % 2 != 0:
                            L[ip+igr, ig+iEv] = L[ip+igr, ig+iEv]+(Neq+1)/(2*Neq+1)*A[count+1, jj]
                            L[ip+igr+M-1, ig+iEv+M-1] = L[ip+igr+M-1, ig+iEv+M-1]+(Neq+1)/(2*Neq+1)*A[count+m//2+1, jj]

                    count = count + 1*(moment > 0)

    op.L = L
    return op


def MarkCoeffs(PN, even=False):
    """
    Evaluate Mark boundary condition coefficients.

    Parameters
    ----------
    PN : int
        Spherical harmonics approximation order.
    even : bool, optional
        If true, previous order Legendre polynomials roots are used instead
        of current even order roots. The default is ``False``.

    Returns
    -------
    A : numpy.ndarray
        Matrix collecting Mark coefficients.

    """
    # FIXME: to correctly handle one-side BC, only positive/negative directions
    # have to be kept!
    if PN % 2 == 0:
        M = PN-1
    else:
        M = PN

    if even is True and PN % 2 == 0:
        roots, weights = roots_legendre(PN)
    else:
        roots, weights = roots_legendre(PN+1)

    roots = np.delete(roots, np.where(roots == 0))
    roots = -np.sort(-roots)  # sort in descending order

    A = np.zeros((M+1, PN+1))
    for row in range(0, M+1):
        normcoeffs = (2*(np.arange(0, PN+1))+1)/2
        A[row, :] = normcoeffs*eval_legendre(np.arange(0, PN+1), roots[row])

    return A


def MarshakCoeffs(PN):
    """
    Evaluate Marshak boundary condition coefficients.

    Parameters
    ----------
    PN : int
        Spherical harmonics approximation order.

    Returns
    -------
    A : numpy.ndarray
        Matrix collecting Mark coefficients.

    """
    if PN % 2 == 0:
        M = PN-1

    else:
        M = PN

    normcoeffs = (2*(np.arange(0, PN+1))+1)/2

    m = 1
    A = np.zeros((M+1, PN+1))

    try:
        path = os.path.join(os.path.dirname(__file__), 'Marshak.txt')
        coeffs = np.loadtxt(path)
        # r, c, coeffs = np.loadtxt(path, unpack=True,
        #                           dtype={'names':
        #                                  ('row', 'col', 'vals'),
        #                                  'formats': (np.int, np.int, np.float)})
        rows, cols = coeffs.shape

        if rows >= (M+1)//2 and cols >= (PN+1)//2:
            A[0:(M+1)//2, :] = coeffs[0:(M+1)//2, 0:PN+1]
            A[(M+1)//2:, 0::2] = -coeffs[0:(M+1)//2, 0:PN+1:2]
            A[(M+1)//2:, 1::2] = coeffs[0:(M+1)//2, 1:PN+1:2]

        else:
            raise OSError

    except OSError:

        x = Symbol('x')
        for row in range(0, (M+1)//2):

            for n in range(0, PN+1):
                # positive range
                if n % 2 == 0:  # even moments
                    # posI, err = quad(lambda x:
                    #                  lpmv(0, n, x)*lpmv(0, m, x), 0, 1)
                    f = legendre(n, x)*legendre(m, x)
                    posI = integrate(f, (x, 0, 1))
                    posI = normcoeffs[n]*posI
                    # negative half range
                    A[row+(M+1)//2, n] = -posI

                else:  # odd moments

                    if n == m:
                        posI = 1/2
                        # negative half range
                        A[row+(M+1)//2, n] = posI

                    else:
                        posI = 0
                        # negative half range
                        A[row+(M+1)//2, n] = 0

                # positive half range
                A[row, n] = posI

            # +2 to get odd moments
            m = m+2

        path = os.path.join(os.path.dirname(__file__), 'Marshak.txt')
        np.savetxt(path, A[0:(M+1)//2, :], '%-.10e', delimiter='  ')

    return A


def _getcoeffs(A):
    """
    Compute boundary condition coefficients.

    Odd moments coefficients are evaluated as a function
    of even moments coefficients. The output B matrix
    columns are thus the even moments coefficients in
    increasing order while its rows are the equations
    associated to each boundary condition (Mark or Marshak).

    Parameters
    ----------
    A : numpy.ndarray
        Matrix collecting boundary condition coefficients.

    Returns
    -------
    B : numpy.ndarray
        Its columns are the even moments coefficients in
        increasing order while its rows are the equations
        associated to each boundary condition equation
        (Mark or Marshak).

    """
    m, n = A.shape
    # positive coefficients
    Apos_even = -A[0:m//2, 0::2]
    Apos_odd = A[0:m//2, 1::2]
    # negative coefficients
    Aneg_even = -A[m//2:m, 0::2]
    Aneg_odd = A[m//2:m, 1::2]

    B = np.concatenate([np.dot(np.linalg.inv(Apos_odd), Apos_even),
                        np.dot(np.linalg.inv(Aneg_odd), Aneg_even)])

    return B
