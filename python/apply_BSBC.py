"""
apply_BSBC.py

Derives the BSBC estimates of the mean responses of growth to natural
disasters in the model of:
  Thomas B. Fomby, Yuki Ikeda, and Norman V. Loayza (2013),
  "The Growth Aftermath of Natural Disasters", Journal of Applied Econometrics.

Translated from MATLAB code by Yuki Ikeda (8/10/2011).
"""

import time
from pathlib import Path
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
import pandas as pd

from helpers_numba import psimat, meanpsi, btst

t_start = time.time()

# Set parameter values
B = 1000                 # number of bootstrap samples
st = 6                   # number of periods for the forecast
dep = 'Mean Response of GDP growth to'

p = 2;  q = 2            # Panel VARX(p,q)

# Load data (dev_grow, dev_agr, dev_non, rich_grow, rich_agr, rich_non)
DATASET = 'dev_grow'
DATA_DIR = Path(__file__).resolve().parent.parent / 'original_package' / 'data'
vars_ = pd.read_csv(DATA_DIR / f'{DATASET}.csv', header=None).values
nobs, col = vars_.shape

# Country id and time dimensions
newid = vars_[:, col - 1].astype(int)
N = int(newid.max())
Ti = np.array([np.sum(newid == i) for i in range(1, N + 1)])

acTi = np.zeros(N + 1, dtype=int)
for i in range(N):
    acTi[i + 1] = acTi[i] + Ti[i]

# Set model configuration
K = 6       # number of endogenous variables
dis = 6     # number of exogenous variables

Y = np.column_stack([
    vars_[:, 144], vars_[:, 93], vars_[:, 145], vars_[:, 146],
    vars_[:, 147], vars_[:, 1]
])

LYs = np.column_stack([
    vars_[:, 148], vars_[:, 103], vars_[:, 149], vars_[:, 150],
    vars_[:, 151], vars_[:, 3],
    vars_[:, 152], vars_[:, 113], vars_[:, 153], vars_[:, 154],
    vars_[:, 155], vars_[:, 5]
])

nv = K * p + dis * (q + 1)   # number of explanatory variables per equation

# Benchmark model [Table III in the paper]
Xs = np.column_stack([
    vars_[:, 8], vars_[:, 9], vars_[:, 10], vars_[:, 11],
    vars_[:, 0], vars_[:, 97],
    vars_[:, 16], vars_[:, 17], vars_[:, 18], vars_[:, 19],
    vars_[:, 2], vars_[:, 107],
    vars_[:, 24], vars_[:, 25], vars_[:, 26], vars_[:, 27],
    vars_[:, 4], vars_[:, 117]
])

# Right-hand-side matrix
RH = np.hstack([LYs, Xs])

# Matrix D (fixed-effect dummies)
D = np.zeros((nobs, N))
for nn in range(N):
    D[acTi[nn]:acTi[nn + 1], nn] = 1.0

# Matrix A (within transformation)
A = np.zeros((nobs, nobs))
for i in range(N):
    s = acTi[i]
    e = acTi[i + 1]
    t = Ti[i]
    A[s:e, s:e] = np.eye(t) - (1.0 / t) * np.ones((t, 1)) @ np.ones((1, t))

# LSDV estimator
deltahat = inv(RH.T @ A @ RH) @ RH.T @ A @ Y

# Bootstrap procedure
# Derive the mean of psi_b
deltatil_b = btst(B, nobs, nv, p, K, N, A, Y, LYs, Xs, Ti, acTi, D, deltahat)
m_psi_b = meanpsi(p, q, K, dis, st, K, B, deltatil_b)

# Calculate psi of Pesaran and Zhao
psi_pz = 2 * psimat(p, q, K, dis, st, deltahat) - m_psi_b

# Store mean responses
mr = np.zeros((st + 1, 4))
for i in range(st + 1):
    mr[i, :] = psi_pz[K * (i + 1) - 1, 0:4]

# Derive BSBC estimate for each bootstrap replicate
mr_store = np.zeros((st + 1, 4 * B))

for h in range(B):
    beta_h = deltatil_b[:, h * K : (h + 1) * K]
    deltatil_b_h = btst(B, nobs, nv, p, K, N, A, Y, LYs, Xs, Ti, acTi, D, beta_h)
    m_psi_b_h = meanpsi(p, q, K, dis, st, K, B, deltatil_b_h)
    psi_pz_h = 2 * psimat(p, q, K, dis, st, beta_h) - m_psi_b_h

    for i in range(st + 1):
        mr_store[i, h * 4 : (h + 1) * 4] = psi_pz_h[K * (i + 1) - 1, 0:4]

    print(h + 1)

# Reshape and sort bootstrap estimates
mr_store_long = np.zeros((B, 4 * (st + 1)))
for h in range(B):
    mr_store_long[h, :] = mr_store[:, h * 4 : (h + 1) * 4].T.ravel()

mr_store_long = np.sort(mr_store_long, axis=0)

mr_low5  = mr_store_long[round(B * 0.05) - 1, :].reshape(4, st + 1).T
mr_low10 = mr_store_long[round(B * 0.1) - 1, :].reshape(4, st + 1).T
mr_up90  = mr_store_long[round(B * 0.90) - 1, :].reshape(4, st + 1).T
mr_up95  = mr_store_long[round(B * 0.95) - 1, :].reshape(4, st + 1).T

# Store cumulative effects
ce_store = mr_store[0:4, :].sum(axis=0)
ce_store_long = np.zeros((B, 4))
for h in range(B):
    ce_store_long[h, :] = ce_store[h * 4 : (h + 1) * 4]
ce_store_long = np.sort(ce_store_long, axis=0)

ce = np.zeros((5, 4))
ce[0, :] = ce_store_long[round(B * 0.05) - 1, :]
ce[1, :] = ce_store_long[round(B * 0.1) - 1, :]
ce[2, :] = mr[0:4, :].sum(axis=0)
ce[3, :] = ce_store_long[round(B * 0.9) - 1, :]
ce[4, :] = ce_store_long[round(B * 0.95) - 1, :]

# Significance helpers
def _sig_mr(val, up90, low10, up95, low5):
    if up90 > 0 and low10 > 0 and low5 < 0:
        return '*'
    elif up95 > 0 and low5 > 0:
        return '**'
    elif up90 < 0 and low10 < 0 and up95 > 0:
        return '*'
    elif up95 < 0 and low5 < 0:
        return '**'
    return ''

def _sig_ce(val, ce_l5, ce_l10, ce_u90, ce_u95):
    if ce_u90 > 0 and ce_l10 > 0 and ce_l5 < 0:
        return '*'
    elif ce_u95 > 0 and ce_l5 > 0:
        return '**'
    elif ce_u90 < 0 and ce_l10 < 0 and ce_u95 > 0:
        return '*'
    elif ce_u95 < 0 and ce_l5 < 0:
        return '**'
    return ''

def _fmt(val, sig):
    s = f'{val:>9.4f}'
    return s + sig

# Print table of mean responses (BSBC estimator)
disaster_names = ['Droughts', 'Earthquakes', 'Floods', 'Storms']

print('\n' + '='*72)
print(f'  Table: {dep} natural disaster shocks')
print('         (BSBC estimator)')
print('='*72)
print(f'{"Disaster":<16} {"Period":<20} {"Estimate":>14}')
print('-'*72)

for j in range(4):
    name = disaster_names[j]
    for i in range(4):
        sig = _sig_mr(mr[i,j], mr_up90[i,j], mr_low10[i,j], mr_up95[i,j], mr_low5[i,j])
        label = name if i == 0 else ''
        print(f'{label:<16} {"Year "+str(i):<20} {_fmt(mr[i,j], sig):>14}')
    sig_c = _sig_ce(ce[2,j], ce[0,j], ce[1,j], ce[3,j], ce[4,j])
    print(f'{"":<16} {"Cumulative effect":<20} {_fmt(ce[2,j], sig_c):>14}')

print('-'*72)
print(f'{"":>36} {"N countries":>14}  {N}')
print(f'{"":>36} {"Sample period":>14}  1960-2007')
print(f'{"":>36} {"Observations":>14}  {nobs}')
print('-'*72)
print('Note: * significant at one-tail 10% level, ** at 5% level.')
print('='*72)

elapsed = time.time() - t_start
print(f"\nElapsed time: {elapsed:.2f} s  ({elapsed / 60:.2f} min)")

# Draw graphs
x = np.arange(st + 1)
y = np.zeros(st + 1)
titles = ['Drought Shock', 'Earthquake Shock', 'Flood Shock', 'Storm Shock']

for j in range(4):
    plt.figure()
    plt.plot(x, mr[:, j], '-b', label='BSBC estimate')
    plt.plot(x, mr_up90[:, j], ':b', label='90% upper')
    plt.plot(x, mr_low10[:, j], ':b', label='10% lower')
    plt.plot(x, y, '-k')
    plt.xlabel('Time (years)')
    plt.title(f'{dep} {titles[j]}')
    plt.legend()
    plt.tight_layout()

plt.show()
