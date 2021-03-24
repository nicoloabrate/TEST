"""
Author: N. Abrate.

File: GPT.py

Description: Class for Generalised Perturbation Theory method.
"""
import numpy as np


class GPT():

    def __init__(self, geom, unperturbed, perturbed, psi, phi, mu, N, M):
        A, B = unperturbed.A, unperturbed.B
        Ap, Bp = perturbed.A, perturbed.B

        # normalisation constants
        C = np.zeros(phi.shape[1],1)
        for i in range(0, M):
            C[i] = GPT.braket(psi[:,i], np.dot(B, psi[:,i]), geom)

        dB = Bp-B
        dA = (Ap-A)-mu[0]*dB
        lambdas, p = GPT.computeperturbations(N, dA, dB, B, Bp, mu, phi, psi, C, geom)

    def computeperturbations(N, dA, dB, B, Bp, mu, fmodes, amodes, C, geom):
        """
        Compute the eigenvalue problem perturbations

        Parameters
        ----------
        N : int
            DESCRIPTION.
        dA : ndarray
            DESCRIPTION.
        dB : ndarray
            DESCRIPTION.
        B : ndarray
            DESCRIPTION.
        Bp : ndarray
            DESCRIPTION.
        mu : ndarray
            DESCRIPTION.
        fmodes : ndarray
            DESCRIPTION.
        amodes : ndarray
            DESCRIPTION.
        C : ndarray
            DESCRIPTION.
        geom : object
            DESCRIPTION.

        Raises
        ------
        OSError
            Forward and adjoint modes dimension mismatch!

        Returns
        -------
        None.

        """
        if fmodes.shape[1] != amodes.shape[1]:
            raise OSError('Forward and adjoint modes dimension mismatch!')
        else:
            M = fmodes.shape[1]
        # pre-allocation
        alpha = np.zeros((M,N-1))  # alpha_{1,K} = 0 as normalization constant
        lambdas = np.zeros((N-1,1))  # perturbed eigenvalue
        p = np.zeros((M, N-1)) # perturbed system flux perturbations
        # normalisation constants
        Pp = GPT.braket(1,np.dot(Bp, p[:, 0]), geom)
        
        # eigenvalue perturbations
        for n in range(0, N):
            
            S = np.dot(dA, fmodes[:, 0]) if n == 0 else np.dot(dA, p[:, n-1])
            S1, S2 = 0, 0
            # sum previous contributions
            for k in range(0, n-1):
                if n > k:
                    S1 = S1+lambdas[k]*np.dot(B, p[:, n-k])
                    if n > k+1:
                        S2 = S2+lambdas[k]*np.dot(dB, p[:, n-k-1])
                    else:
                        S2 = S2+lambdas[k]*np.dot(dB, fmodes[:, 0])

            lambdas[n] = (GPT.braket(amodes[:, 0], S, geom)+
                        -GPT.braket(amodes[:, 0], S1, geom)+
                        -GPT.braket(amodes[:, 0], S2, geom))/C[0]
                        
            # eigenvector perturbations
            for m in range(1, M):
                # initialisation 
                S1a, S2a = 0, 0
                for k in range(0, n-1):
                    if n > k:
                        S1a = S1a+lambdas[k]*alpha[m, n-k]
                        
                        if n > k+1:
                            S2a = S2a+lambdas[k]*dB*p[:, n-k-1]
                        else:
                            S2a = S2a+lambdas[k]*dB*fmodes[:, 0]
            
                # coefficient evaluation
                alpha[m, n] = 1/((mu[m]-mu[0])*C[m])*(-GPT.braket(amodes[:, m], S, geom)+
                                                      +S1a*C[m]+GPT.braket(amodes[:, m], S2a, geom))
                
                if n == 0:
                    alpha[0, n] = (-GPT.braket(1, np.dot(dB, fmodes[:, 0]), geom)+GPT.braket(1, Bp*np.dot(fmodes[:, 1:M], alpha[1:M, n]), geom))/Pp
                else:
                    alpha[0, n] = -GPT.braket(1, Bp*np.dot(fmodes[:, 1:M], alpha[1:M, n]), geom)/Pp
                # perturbation reconstruction
                p[:, n] = np.dot(fmodes[:, 0:M], alpha[:, n])
         
        return lambdas, p

    def braket(v1, v2, geom):
        """
        Compute bra-ket product over the phase space.

        Parameters
        ----------
        v1 : ndarray
            DESCRIPTION.
        v2 : ndarray
            DESCRIPTION.
        geom : object
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if geom.nA > 0:
            raise OSError('GPT cannot be applied yet to transport solutions!')

        v1v2 = np.multiply(v1, v2)
        G = geom.nE
        S = geom.nS
        grid = geom.mesh
        I = 0
        for g in range(0, G):
            skip = g*S
            I = np.trapz(v1v2[skip:skip+S], x=grid)
        return I
    
            