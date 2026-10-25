"""
helpers_numba.py

Numba JIT-compiled versions of psimat, meanpsi, and btst for Panel VARX estimation.
"""

import numpy as np
from numba import njit


@njit(cache=True)
def psimat(p, q, nen, dis, st, beta):
    """Calculate impulse-response matrix psi from phi and theta."""
    # Set up phi matrix
    phi = np.zeros((nen * p, nen))
    for i in range(p):
        phi[nen * i : nen * (i + 1), :] = beta[nen * i : nen * (i + 1), :].T

    phi1 = phi[0:nen, :].copy()
    if p >= 2:
        phi2 = phi[nen : 2 * nen, :].copy()
    else:
        phi2 = np.zeros((nen, nen))
    phi3 = np.zeros((nen, nen))

    # Set up theta matrix
    theta = np.zeros((nen * (q + 1), dis))
    for i in range(q + 1):
        theta[i * nen : (i + 1) * nen, :] = beta[nen * p + dis * i : nen * p + dis * (i + 1), :].T

    theta0 = theta[0:nen, :].copy()
    theta1 = theta[nen : 2 * nen, :].copy()
    theta2 = theta[2 * nen : 3 * nen, :].copy()
    theta3 = np.zeros((nen, dis))

    # Calculate psi matrix
    psi = np.zeros((nen * (st + 1), dis))
    psi[0:nen, :] = theta0
    psi[nen : 2 * nen, :] = theta1 + phi1 @ np.ascontiguousarray(psi[0:nen, :])
    psi[2 * nen : 3 * nen, :] = (theta2
        + phi1 @ np.ascontiguousarray(psi[nen : 2 * nen, :])
        + phi2 @ np.ascontiguousarray(psi[0:nen, :]))
    psi[3 * nen : 4 * nen, :] = (theta3
        + phi1 @ np.ascontiguousarray(psi[2 * nen : 3 * nen, :])
        + phi2 @ np.ascontiguousarray(psi[nen : 2 * nen, :])
        + phi3 @ np.ascontiguousarray(psi[0:nen, :]))

    for i in range(4, st + 1):
        psi[nen * i : nen * (i + 1), :] = (
            phi1 @ np.ascontiguousarray(psi[nen * (i - 1) : nen * i, :])
            + phi2 @ np.ascontiguousarray(psi[nen * (i - 2) : nen * (i - 1), :])
            + phi3 @ np.ascontiguousarray(psi[nen * (i - 3) : nen * (i - 2), :])
        )

    return psi


@njit(cache=True)
def meanpsi(p, q, nen, dis, st, K, B, deltatil_b):
    """Average of psi_b over B bootstrap replications."""
    m_psi_b = np.zeros((nen * (st + 1), dis))

    for i in range(B):
        m_psi_b += psimat(p, q, nen, dis, st,
                          deltatil_b[:, i * K : (i + 1) * K])

    m_psi_b /= B
    return m_psi_b


@njit(cache=True)
def btst(B, nobs, nv, p, K, N, A, Y, LYs, Xs, Ti, acTi, D, deltahat):
    """Bootstrap procedure for calculating deltatil_b."""
    RH = np.hstack((LYs, Xs))
    deltatil_b = np.zeros((nv, K * B))
    deltatil_old = deltahat.copy()

    # Calculate alphatil
    residual = Y - RH @ deltatil_old
    alphatil_pre = D.T @ residual
    alphatil = np.zeros((N, K))
    for i in range(N):
        alphatil[i, :] = alphatil_pre[i, :] / Ti[i]
    D_alphatil = D @ alphatil

    # Residuals
    epstil = Y - RH @ deltatil_old - D @ alphatil

    # Rescaling residuals
    ARH = A @ RH
    RtARH = RH.T @ ARH
    RtA = RH.T @ A
    X_solve = np.linalg.solve(RtARH, RtA)
    proj = ARH @ X_solve

    M_diag = np.zeros(nobs)
    for i in range(nobs):
        M_diag[i] = A[i, i] - proj[i, i]

    eom = np.zeros((nobs, K))
    for i in range(nobs):
        sqrt_m = np.sqrt(M_diag[i]) if M_diag[i] > 0 else 1.0
        for k in range(K):
            eom[i, k] = epstil[i, k] / sqrt_m

    epstil_r = np.zeros((nobs, K))
    for i in range(N):
        s = acTi[i]
        e = acTi[i + 1]
        t = Ti[i]
        block_mean = np.zeros(K)
        for row in range(s, e):
            for k in range(K):
                block_mean[k] += eom[row, k]
        for k in range(K):
            block_mean[k] /= t
        scale = np.sqrt(float(t) / float(t - 1))
        for row in range(s, e):
            for k in range(K):
                epstil_r[row, k] = scale * (eom[row, k] - block_mean[k])

    nXs = Xs.shape[1]

    # Bootstrap sampling (B times)
    for brep in range(B):
        rb = np.random.randint(0, nobs, size=nobs)
        epstil_b = epstil_r[rb, :]

        # Set initial lagged values
        ly_b = np.zeros((nobs, p * K))
        for i in range(N):
            s = acTi[i]
            for k in range(p * K):
                ly_b[s, k] = LYs[s, k]
            if p >= 2 and s + 1 < acTi[i + 1]:
                for k in range(K, 2 * K):
                    ly_b[s + 1, k] = LYs[s + 1, k]

        # Generate y_b recursively
        y_b = np.zeros((nobs, K))
        for nn in range(N):
            s = acTi[nn]
            e = acTi[nn + 1]
            for i in range(s, e):
                for k in range(K):
                    val = D_alphatil[i, k] + epstil_b[i, k]
                    for j in range(p * K):
                        val += ly_b[i, j] * deltatil_old[j, k]
                    for j in range(nXs):
                        val += Xs[i, j] * deltatil_old[p * K + j, k]
                    y_b[i, k] = val

                # Update lagged values
                if i <= e - 3:
                    for k in range(K):
                        ly_b[i + 1, k] = y_b[i, k]
                        ly_b[i + 2, K + k] = y_b[i, k]
                elif i == e - 2:
                    for k in range(K):
                        ly_b[i + 1, k] = y_b[i, k]

        # LSDV estimator for bootstrap sample
        W_b = np.hstack((ly_b, Xs))
        WtA = W_b.T @ A
        WtAW = WtA @ W_b
        WtAy = WtA @ y_b
        result = np.linalg.solve(WtAW, WtAy)
        deltatil_b[:, brep * K : (brep + 1) * K] = result

    return deltatil_b
