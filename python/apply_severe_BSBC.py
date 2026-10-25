"""
apply_severe_BSBC.py

Derives the BSBC estimates of the mean responses of growth to natural
disasters, with severe disasters differentiated.

Reference:
  Fomby, Ikeda, and Loayza (2013), "The Growth Aftermath of Natural Disasters", JAE.

Translated from MATLAB code by Yuki Ikeda (8/10/2011).
"""

import time, numpy as np
from pathlib import Path
from numpy.linalg import inv
import matplotlib.pyplot as plt, pandas as pd
from helpers_numba import psimat, meanpsi, btst

t_start = time.time()

# Parameters
B = 1000; st = 6; dep = 'Mean Response of GDP growth to (severe)'
p = 2; q = 2

# Load data (dev_grow, dev_agr, dev_non)
DATASET = 'dev_grow'
DATA_DIR = Path(__file__).resolve().parent.parent / 'original_package' / 'data'
vars_ = pd.read_csv(DATA_DIR / f'{DATASET}.csv', header=None).values
nobs, col = vars_.shape
newid = vars_[:, col-1].astype(int); N = int(newid.max())
Ti = np.array([np.sum(newid == i) for i in range(1, N+1)])
acTi = np.zeros(N+1, dtype=int)
for i in range(N): acTi[i+1] = acTi[i] + Ti[i]

# Model configuration
K = 6; dis = 10
Y = np.column_stack([vars_[:,144], vars_[:,93], vars_[:,145], vars_[:,146], vars_[:,147], vars_[:,1]])
LYs = np.column_stack([vars_[:,148], vars_[:,103], vars_[:,149], vars_[:,150], vars_[:,151], vars_[:,3],
                        vars_[:,152], vars_[:,113], vars_[:,153], vars_[:,154], vars_[:,155], vars_[:,5]])
# Benchmark model [Table V in the paper]
Xs = np.column_stack([vars_[:,8], vars_[:,9], vars_[:,10], vars_[:,11],
                       vars_[:,12], vars_[:,13], vars_[:,14], vars_[:,15],
                       vars_[:,0], vars_[:,97],
                       vars_[:,16], vars_[:,17], vars_[:,18], vars_[:,19],
                       vars_[:,20], vars_[:,21], vars_[:,22], vars_[:,23],
                       vars_[:,2], vars_[:,107],
                       vars_[:,24], vars_[:,25], vars_[:,26], vars_[:,27],
                       vars_[:,28], vars_[:,29], vars_[:,30], vars_[:,31],
                       vars_[:,4], vars_[:,117]])

nv = K*p + dis*(q+1)
RH = np.hstack([LYs, Xs])

# Matrices D and A
D = np.zeros((nobs,N))
for nn in range(N): D[acTi[nn]:acTi[nn+1], nn] = 1.0
A = np.zeros((nobs,nobs))
for i in range(N):
    s,e,t = acTi[i], acTi[i+1], Ti[i]
    A[s:e,s:e] = np.eye(t) - np.ones((t,t))/t

# LSDV estimator
deltahat = inv(RH.T @ A @ RH) @ RH.T @ A @ Y

# Bootstrap procedure
deltatil_b = btst(B, nobs, nv, p, K, N, A, Y, LYs, Xs, Ti, acTi, D, deltahat)
m_psi_b = meanpsi(p, q, K, dis, st, K, B, deltatil_b)
psi_pz = 2*psimat(p,q,K,dis,st,deltahat) - m_psi_b

# Mean responses
mr = np.zeros((st+1, 8))
for i in range(st+1): mr[i,:] = psi_pz[K*(i+1)-1, 0:8]

mr_sev = np.zeros((st+1, 4))
for i in range(4): mr_sev[:,i] = mr[:,i] + mr[:,i+4]

# BSBC for each bootstrap replicate
mr_store = np.zeros((st+1, 8*B)); mr_sev_store = np.zeros((st+1, 4*B))
for h in range(B):
    beta_h = deltatil_b[:, h*K:(h+1)*K]
    deltatil_b_h = btst(B, nobs, nv, p, K, N, A, Y, LYs, Xs, Ti, acTi, D, beta_h)
    m_psi_b_h = meanpsi(p,q,K,dis,st,K,B,deltatil_b_h)
    psi_pz_h = 2*psimat(p,q,K,dis,st,beta_h) - m_psi_b_h
    for i in range(st+1):
        mr_store[i, h*8:(h+1)*8] = psi_pz_h[K*(i+1)-1, 0:8]
        mr_sev_store[i, h*4:(h+1)*4] = mr_store[i, h*8:h*8+4] + mr_store[i, h*8+4:(h+1)*8]
    print(h+1)

# Reshape and sort
mr_store_long = np.zeros((B, 8*(st+1))); mr_sev_store_long = np.zeros((B, 4*(st+1)))
for h in range(B):
    mr_store_long[h,:] = mr_store[:, h*8:(h+1)*8].T.ravel()
    mr_sev_store_long[h,:] = mr_sev_store[:, h*4:(h+1)*4].T.ravel()

mr_store_long = np.sort(mr_store_long, axis=0)
mr_low5  = mr_store_long[round(B*0.05)-1,:].reshape(8,st+1).T
mr_low10 = mr_store_long[round(B*0.1)-1,:].reshape(8,st+1).T
mr_up90  = mr_store_long[round(B*0.90)-1,:].reshape(8,st+1).T
mr_up95  = mr_store_long[round(B*0.95)-1,:].reshape(8,st+1).T

mr_sev_store_long = np.sort(mr_sev_store_long, axis=0)
mr_sev_low5  = mr_sev_store_long[round(B*0.05)-1,:].reshape(4,st+1).T
mr_sev_low10 = mr_sev_store_long[round(B*0.1)-1,:].reshape(4,st+1).T
mr_sev_up90  = mr_sev_store_long[round(B*0.90)-1,:].reshape(4,st+1).T
mr_sev_up95  = mr_sev_store_long[round(B*0.95)-1,:].reshape(4,st+1).T

# Cumulative effects
ce_store = mr_store[0:4,:].sum(axis=0); ce_store_long = np.zeros((B,8))
ce_sev_store = mr_sev_store[0:4,:].sum(axis=0); ce_sev_store_long = np.zeros((B,4))
for h in range(B):
    ce_store_long[h,:] = ce_store[h*8:(h+1)*8]
    ce_sev_store_long[h,:] = ce_sev_store[h*4:(h+1)*4]
ce_store_long = np.sort(ce_store_long, axis=0)
ce_sev_store_long = np.sort(ce_sev_store_long, axis=0)

ce = np.zeros((5,8))
ce[0] = ce_store_long[round(B*0.05)-1]; ce[1] = ce_store_long[round(B*0.1)-1]
ce[2] = mr[0:4,:].sum(0); ce[3] = ce_store_long[round(B*0.9)-1]; ce[4] = ce_store_long[round(B*0.95)-1]

ce_sev = np.zeros((5,4))
ce_sev[0] = ce_sev_store_long[round(B*0.05)-1]; ce_sev[1] = ce_sev_store_long[round(B*0.1)-1]
ce_sev[2] = mr_sev[0:4,:].sum(0); ce_sev[3] = ce_sev_store_long[round(B*0.9)-1]
ce_sev[4] = ce_sev_store_long[round(B*0.95)-1]

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

disaster_names = ['Droughts', 'Earthquakes', 'Floods', 'Storms']

# Combined table: Base vs Severe disasters (BSBC estimator)
print('\n' + '='*82)
print(f'  Table: {dep} natural disaster shocks')
print('         (BSBC estimator - base vs. severe)')
print('='*82)
print(f'{"Disaster":<16} {"Period":<20} {"Base":>14} {"Severe":>14}')
print('-'*82)

for j in range(4):
    name = disaster_names[j]
    for i in range(4):
        sig_b = _sig_mr(mr[i,j], mr_up90[i,j], mr_low10[i,j], mr_up95[i,j], mr_low5[i,j])
        sig_s = _sig_mr(mr_sev[i,j], mr_sev_up90[i,j], mr_sev_low10[i,j],
                         mr_sev_up95[i,j], mr_sev_low5[i,j])
        label = name if i == 0 else ''
        print(f'{label:<16} {"Year "+str(i):<20} {_fmt(mr[i,j], sig_b):>14} {_fmt(mr_sev[i,j], sig_s):>14}')
    sig_cb = _sig_ce(ce[2,j], ce[0,j], ce[1,j], ce[3,j], ce[4,j])
    sig_cs = _sig_ce(ce_sev[2,j], ce_sev[0,j], ce_sev[1,j], ce_sev[3,j], ce_sev[4,j])
    print(f'{"":<16} {"Cumulative effect":<20} {_fmt(ce[2,j], sig_cb):>14} {_fmt(ce_sev[2,j], sig_cs):>14}')

print('-'*82)
print(f'{"":>36} {"N countries":>14}  {N}')
print(f'{"":>36} {"Sample period":>14}  1960-2007')
print(f'{"":>36} {"Observations":>14}  {nobs}')
print('-'*82)
print('Note: * significant at one-tail 10% level, ** at 5% level.')
print('      Base = moderate disasters; Severe = moderate + severe (total effect).')
print('='*82)

elapsed = time.time() - t_start
print(f"\nElapsed: {elapsed:.1f}s ({elapsed/60:.2f}min)")

# Graphs
x = np.arange(st+1); y = np.zeros(st+1)
for j, t in enumerate(['Drought','Earthquake','Flood','Storm']):
    plt.figure()
    plt.plot(x, mr[:,j], '-b', linewidth=2, label='Base')
    plt.plot(x, mr_up90[:,j], ':b', alpha=0.5)
    plt.plot(x, mr_low10[:,j], ':b', alpha=0.5)
    plt.plot(x, mr_sev[:,j], '-r', linewidth=2, label='Severe')
    plt.plot(x, mr_sev_up90[:,j], ':r', alpha=0.5)
    plt.plot(x, mr_sev_low10[:,j], ':r', alpha=0.5)
    plt.plot(x, y, '-k', linewidth=0.5)
    plt.xlabel('Time (years)')
    plt.title(f'{dep} {t} Shock')
    plt.legend()
    plt.tight_layout()
plt.show()
