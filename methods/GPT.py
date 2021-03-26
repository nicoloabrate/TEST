"""
Author: N. Abrate.

File: GPT.py

Description: Class for Generalised Perturbation Theory method.
"""
from numpy import zeros, dot
from TEST.phasespace import PhaseSpace
from TEST.phasespace import interp


class GPT():

    def __init__(self, N, geom, unperturbed, pertgeom, perturbed, fwd, adj):

        self.PertOrder = N
        # check forward and adjoint modes consistency
        M, M1 = len(fwd.eigvals), len(adj.eigvals)
        if M != M1:
            raise OSError('Forward (%d) and adjoint (%d) modes mismatch!' % (M, M1))
        else:
            self.ModExpOrder = M
            phi, psi =  fwd.eigvect, adj.eigvect
            mu = fwd.eigvals

        # define phase space
        PS = PhaseSpace(geom)
        xp = pertgeom.mesh

        # --- get operators
        A, B = unperturbed.A, unperturbed.B
        Ap, Bp = perturbed.A, perturbed.B

        # normalisation constants
        C = zeros(phi.shape[1],1)
        for i in range(0, self.ModExpOrder):
            C[i] = PS.braket(psi[:,i], dot(B, psi[:,i]))

        # pre-allocation
        alpha = zeros((M, N-1))  # alpha_{1,K} = 0 as normalization constant
        lambdas = zeros((N-1, 1))  # perturbed eigenvalue
        p = zeros((M, N-1)) # perturbed system flux perturbations
        # normalisation constants
        tmp = PS.interp(p[:, 0], xp)  # interpolate over pert. mesh
        tmp = PS.interp(dot(Bp, tmp))  # interpolate product over ref. mesh
        Pp = PS.braket(tmp)

        # eigenvalue perturbations
        for n in range(0, self.N):
            if n == 0:
                v1 = PS.interp(phi[:, 0], xp)
                S = PS.interp(dot(Ap-mu[0]*Bp, v1))-dot(A-mu[0]*B, phi[:, 0])
            else:
                v1 = PS.interp(p[:, n-1], xp)
                S = PS.interp(dot(Ap-mu[0]*Bp, v1))-dot(A-mu[0]*B, p[:, n-1])

            S1, S2 = 0, 0

            # sum previous contributions
            for k in range(0, n-1):
                if n > k:
                    if n > k+1:
                        v1 = PS.interp(p[:, n-k-1], xp)
                        v2 = p[:, n-k-1]
                    else:
                        v1 = PS.interp(phi[:, 0], xp)
                        v2 = phi[:, 0]

                    S1 = S1+lambdas[k]*dot(B, p[:, n-k])
                    S2 = S2+lambdas[k]*(PS.interp(dot(Bp, v1))-dot(B, v2))

            lambdas[n] = (PS.braket(psi[:, 0], S)-PS.braket(psi[:, 0], S1)+
                         -PS.braket(psi[:, 0], S2))/C[0]

            # eigenvector perturbations
            for m in range(1, M):
                # initialisation 
                S1a, S2a = 0, 0
                for k in range(0, n-1):
                    if n > k:
                        if n > k+1:
                            v1 = PS.interp(p[:, n-k-1], xp)
                            v2 = p[:, n-k-1]
                        else:
                            v1 = PS.interp(phi[:, 0], xp)
                            v2 = phi[:, 0]

                        S1a = S1a+lambdas[k]*alpha[m, n-k]
                        S2a = S2a+lambdas[k]*(PS.interp(dot(Bp, v1))-dot(B, v2))

                # coefficient evaluation
                alpha[m, n] = (-PS.braket(psi[:, m], S)+S1a*C[m]+
                               +PS.braket(psi[:, m], S2a))
                alpha[m, n] = alpha[m, n]/((mu[m]-mu[0])*C[m])
                
                if n == 0:
                    v1 = PS.interp(phi[:, 0], xp)
                    Q1 = -PS.braket(PS.interp(dot(Bp, v1))-dot(B, phi[:, 0]))
                    v1 = PS.interp(dot(phi[:, 1:M], alpha[1:M, n]), xp)
                    Q2 = PS.braket(PS.interp(Bp*v1))
                    alpha[0, n] = (Q1+Q2)/Pp
                else:
                    v1 = PS.interp(dot(phi[:, 1:M], alpha[1:M, n]), xp)
                    alpha[0, n] = -PS.braket(PS.interp(Bp*v1))/Pp
                # perturbation reconstruction
                p[:, n] = dot(phi[:, 0:M], alpha[:, n])

        self.EvalPert = lambdas
        self.EvecPert = p
