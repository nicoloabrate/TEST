"""
Author: N. Abrate.

File: GPT.py

Description: Class for Generalised Perturbation Theory method.
"""
from numpy import zeros, dot
from TEST.phasespace import PhaseSpace


class GPT():

    def __init__(self, N, geom, unperturbed, pertgeom, perturbed, adj):

        self.PertOrder = N
        # check forward and adjoint modes consistency
        M, M1 = len(unperturbed.eigvals), len(adj.eigvals)
        if M != M1:
            raise OSError('Forward (%d) and adjoint (%d) modes mismatch!' % (M, M1))
        else:
            self.ModExpOrder = M
            phi, psi =  unperturbed.eigvect, adj.eigvect
            mu = 1/unperturbed.eigvals

        # define phase space
        PS = PhaseSpace(geom)
        xp = pertgeom.mesh

        # --- get operators
        A, B = unperturbed.A, unperturbed.B
        Ap, Bp = perturbed.A, perturbed.B

        # normalisation constants
        C = zeros((phi.shape[1], ))
        for i in range(0, self.ModExpOrder):
            C[i] = PS.braket(psi[:,i], B.dot(psi[:,i]))

        # pre-allocation
        alpha = zeros((M, N))  # alpha_{1,K} = 0 as normalization constant
        lambdas = zeros((N, ))  # perturbed eigenvalue
        p = zeros((phi.shape[0], N)) # perturbed system flux perturbations
        # normalisation constants
        tmp = PS.interp(phi[:, 0], xp, isref=False)  # interpolate over pert. mesh
        tmp = PS.interp(Bp.dot(tmp), xp)  # interpolate product over refer. mesh
        Pp = PS.braket(tmp)

        # eigenvalue perturbations
        for n in range(0, self.PertOrder):
            if n == 0:
                v1 = PS.interp(phi[:, 0], xp, isref=False)
                S = PS.interp((Ap-mu[0]*Bp).dot(v1), xp)-(A-mu[0]*B).dot(phi[:, 0])
            else:
                v1 = PS.interp(p[:, n-1], xp, isref=False)
                S = PS.interp((Ap-mu[0]*Bp).dot(v1), xp)-(A-mu[0]*B).dot(p[:, n-1])

            S1, S2 = 0, 0

            # sum previous contributions
            for k in range(1, n-1):
                if n > k:
                    if n > k+1:
                        v1 = PS.interp(p[:, n-k-1], xp, isref=False)
                        v2 = p[:, n-k-1]
                    else:
                        v1 = PS.interp(phi[:, 0], xp, isref=False)
                        v2 = phi[:, 0]

                    S1 = S1+lambdas[k]*B.dot(p[:, n-k])
                    S2 = S2+lambdas[k]*(PS.interp(Bp.dot(v1), xp)-B.dot(v2))

            lambdas[n] = (PS.braket(psi[:, 0], S)-PS.braket(psi[:, 0], S1)+
                          -PS.braket(psi[:, 0], S2))/C[0]

            # eigenvector perturbations
            for m in range(1, M):
                # initialisation 
                S1a, S2a = 0, 0
                for k in range(1, n-1):
                    if n > k:
                        if n > k+1:
                            v1 = PS.interp(p[:, n-k-1], xp, isref=False)
                            v2 = p[:, n-k-1]
                        else:
                            v1 = PS.interp(phi[:, 0], xp, isref=False)
                            v2 = phi[:, 0]

                        S1a = S1a+lambdas[k]*alpha[m, n-k]
                        S2a = S2a+lambdas[k]*(PS.interp(Bp.dot(v1), xp)-B.dot(v2))

                # coefficient evaluation
                alpha[m, n] = -PS.braket(psi[:, m], S)+S1a*C[m]+PS.braket(psi[:, m], S2a)
                alpha[m, n] = alpha[m, n]/((mu[m]-mu[0])*C[m])
                
                if n == 0:
                    v1 = PS.interp(phi[:, 0], xp, isref=False)
                    Q1 = -PS.braket(PS.interp(Bp.dot(v1), xp)-B.dot(phi[:, 0]))
                    v1 = PS.interp(dot(phi[:, 1:M], alpha[1:M, n]), xp, isref=False)
                    Q2 = PS.braket(PS.interp(Bp.dot(v1), xp))
                    alpha[0, n] = (Q1+Q2)/Pp
                else:
                    v1 = PS.interp(dot(phi[:, 1:M], alpha[1:M, n]), xp, isref=False)
                    alpha[0, n] = -PS.braket(PS.interp(Bp.dot(v1), xp))/Pp
                # perturbation reconstruction
                p[:, n] = dot(phi[:, 0:M], alpha[:, n])

        self.EvalPert = lambdas
        self.ExpCoeff = alpha
        self.EvecPert = p
        self.pertEigv = mu[0]+lambdas.sum()
        self.pertEvec = phi[:, 0]+p.sum(axis=1)
        