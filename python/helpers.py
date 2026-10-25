"""
helpers.py

Helper functions for Panel VARX estimation (psimat, meanpsi, btst).
Translated from MATLAB code by Yuki Ikeda (8/10/2011).
"""

import numpy as np
from numpy.linalg import inv


def psimat(p, q, nen, dis, st, beta):
    """Calculate impulse-response matrix psi from phi and theta."""
    # Set up phi matrix
    phi = np.zeros((nen * p, nen))
    for i in range(p):
        phi[nen * i : nen * (i + 1), :] = beta[nen * i : nen * (i + 1), :].T

    phi1 = phi[0:nen, :]
    phi2 = phi[nen : 2 * nen, :] if p >= 2 else np.zeros((nen, nen))
    phi3 = np.zeros((nen, nen))

    # Set up theta matrix
    theta = np.zeros((nen * (q + 1), dis))
    for i in range(q + 1):
        theta[i * nen : (i + 1) * nen, :] = beta[nen * p + dis * i : nen * p + dis * (i + 1), :].T

    theta0 = theta[0:nen, :]
    theta1 = theta[nen : 2 * nen, :]
    theta2 = theta[2 * nen : 3 * nen, :]
    theta3 = np.zeros((nen, dis))

    # Calculate psi matrix
    psi = np.zeros((nen * (st + 1), dis))
    psi[0:nen, :] = theta0
    psi[nen : 2 * nen, :] = theta1 + phi1 @ psi[0:nen, :]
    psi[2 * nen : 3 * nen, :] = (theta2 + phi1 @ psi[nen : 2 * nen, :]
                                        + phi2 @ psi[0:nen, :])
    psi[3 * nen : 4 * nen, :] = (theta3 + phi1 @ psi[2 * nen : 3 * nen, :]
                                        + phi2 @ psi[nen : 2 * nen, :]
                                        + phi3 @ psi[0:nen, :])

    for i in range(4, st + 1):
        psi[nen * i : nen * (i + 1), :] = (
            phi1 @ psi[nen * (i - 1) : nen * i, :]
            + phi2 @ psi[nen * (i - 2) : nen * (i - 1), :]
            + phi3 @ psi[nen * (i - 3) : nen * (i - 2), :]
        )

    return psi


def meanpsi(p, q, nen, dis, st, K, B, deltatil_b):
    """Average of psi_b over B bootstrap replications."""
    m_psi_b = np.zeros((nen * (st + 1), dis))

    for i in range(B):
        m_psi_b += psimat(p, q, nen, dis, st,
                          deltatil_b[:, i * K : (i + 1) * K])

    m_psi_b /= B
    return m_psi_b


def btst(B, nobs, nv, p, K, N, A, Y, LYs, Xs, Ti, acTi, D, deltahat):
    """Bootstrap procedure for calculating deltatil_b."""
    RH = np.hstack([LYs, Xs])
    deltatil_b = np.zeros((nv, K * B))
    deltatil_old = deltahat.copy()

    # Calculate alphatil
    alphatil_pre = D.T @ (Y - RH @ deltatil_old)
    alphatil = np.zeros((N, K))
    for i in range(N):
        alphatil[i, :] = alphatil_pre[i, :] / Ti[i]
    D_alphatil = D @ alphatil

    # Residuals
    epstil = Y - RH @ deltatil_old - D @ alphatil

    # Rescaling residuals
    ARH = A @ RH
    M_diag = np.diag(A - ARH @ inv(RH.T @ ARH) @ RH.T @ A)
    eom = epstil / np.sqrt(M_diag)[:, None]

    epstil_r = np.zeros((nobs, K))
    for i in range(N):
        s = acTi[i]
        e = acTi[i + 1]
        t = Ti[i]
        block = eom[s:e, :]
        epstil_r[s:e, :] = np.sqrt(t / (t - 1)) * (block - block.mean(axis=0))

    # Bootstrap sampling (B times)
    for brep in range(B):
        rb = np.random.randint(0, nobs, size=nobs)
        epstil_b = epstil_r[rb, :]

        # Set initial lagged values
        ly_b = np.zeros((nobs, p * K))
        for i in range(N):
            s = acTi[i]
            ly_b[s, :] = LYs[s, :]
            if p >= 2:
                ly_b[s + 1, K : 2 * K] = LYs[s + 1, K : 2 * K]

        # Generate y_b recursively
        y_b = np.zeros((nobs, K))
        for nn in range(N):
            s = acTi[nn]
            e = acTi[nn + 1]
            for i in range(s, e):
                rhs = np.hstack([ly_b[i, :], Xs[i, :]])
                y_b[i, :] = rhs @ deltatil_old + D_alphatil[i, :] + epstil_b[i, :]

                # Update lagged values
                if i <= e - 3:
                    ly_b[i + 1, 0:K] = y_b[i, :]
                    ly_b[i + 2, K : 2 * K] = y_b[i, :]
                elif i == e - 2:
                    ly_b[i + 1, 0:K] = y_b[i, :]

        # LSDV estimator for bootstrap sample
        W_b = np.hstack([ly_b, Xs])
        deltatil_b[:, brep * K : (brep + 1) * K] = inv(W_b.T @ A @ W_b) @ W_b.T @ A @ y_b

    return deltatil_b
