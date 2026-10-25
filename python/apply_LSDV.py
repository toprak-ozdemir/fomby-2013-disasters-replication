"""
apply_LSDV.py

Derives the LSDV estimates of the mean responses of growth to natural
disasters using Monte Carlo simulation.

Reference:
  Fomby, Ikeda, and Loayza (2013), "The Growth Aftermath of Natural Disasters", JAE.

Translated from MATLAB code by Yuki Ikeda (8/10/2011).
"""

import time, numpy as np
from pathlib import Path
from numpy.linalg import inv, det
import matplotlib.pyplot as plt, pandas as pd
from helpers import psimat

t_start = time.time()

# Parameters
monte = 30000; st = 6; dep = 'Mean Response of GDP growth to'
p = 2; q = 2

# Load data (dev_grow, dev_agr, dev_non, rich_grow, rich_agr, rich_non)
DATASET = 'dev_grow'
DATA_DIR = Path(__file__).resolve().parent.parent / 'original_package' / 'data'
vars_ = pd.read_csv(DATA_DIR / f'{DATASET}.csv', header=None).values
nobs, col = vars_.shape
newid = vars_[:, col-1].astype(int); N = int(newid.max())
Ti = np.array([np.sum(newid == i) for i in range(1, N+1)])
acTi = np.zeros(N+1, dtype=int)
for i in range(N): acTi[i+1] = acTi[i] + Ti[i]

# Model configuration
K = 6; dis = 6
Y = np.column_stack([vars_[:,144], vars_[:,93], vars_[:,145], vars_[:,146], vars_[:,147], vars_[:,1]])
LYs = np.column_stack([vars_[:,148], vars_[:,103], vars_[:,149], vars_[:,150], vars_[:,151], vars_[:,3],
                        vars_[:,152], vars_[:,113], vars_[:,153], vars_[:,154], vars_[:,155], vars_[:,5]])
nv = K*p + dis*(q+1)
Xs = np.column_stack([vars_[:,8], vars_[:,9], vars_[:,10], vars_[:,11], vars_[:,0], vars_[:,97],
                       vars_[:,16], vars_[:,17], vars_[:,18], vars_[:,19], vars_[:,2], vars_[:,107],
                       vars_[:,24], vars_[:,25], vars_[:,26], vars_[:,27], vars_[:,4], vars_[:,117]])

RH = np.hstack([LYs, Xs])

# Demean
dmY = np.zeros((nobs,K)); dmLYs = np.zeros((nobs,p*K)); dmXs = np.zeros((nobs,(q+1)*dis))
for nn in range(N):
    s, e = acTi[nn], acTi[nn+1]
    dmY[s:e] = Y[s:e] - Y[s:e].mean(0); dmLYs[s:e] = LYs[s:e] - LYs[s:e].mean(0)
    dmXs[s:e] = Xs[s:e] - Xs[s:e].mean(0)
dmRH = np.hstack([dmLYs, dmXs])

# Matrices D and A
D = np.zeros((nobs,N))
for nn in range(N): D[acTi[nn]:acTi[nn+1], nn] = 1.0
A = np.zeros((nobs,nobs))
for i in range(N):
    s,e,t = acTi[i], acTi[i+1], Ti[i]
    A[s:e,s:e] = np.eye(t) - np.ones((t,t))/t

# LSDV estimator
deltahat = inv(RH.T @ A @ RH) @ RH.T @ A @ Y
iQ = inv((dmRH.T @ dmRH)/nobs)
eps = dmY - dmRH @ deltahat
omega = (eps.T @ eps)/nobs
detomega = det(omega)

# Information criteria
like = -(nobs/2.0)*(K*(1+np.log(2*np.pi))+np.log(detomega))
AIC = -2*(like/nobs) + 2*(K*nv/nobs)
SBC = -2*(like/nobs) + (K*nv)*np.log(nobs)/nobs

# Psi from deltahat
psi = psimat(p,q,K,dis,st,deltahat)

# Phi and theta
phi = np.zeros((K*p,K))
for i in range(p): phi[K*i:K*(i+1),:] = deltahat[K*i:K*(i+1),:].T
theta = np.zeros((K*(q+1),dis))
for i in range(q+1): theta[i*K:(i+1)*K,:] = deltahat[K*p+dis*i:K*p+dis*(i+1),:].T

mr = np.zeros((st+1,dis))
for i in range(st+1): mr[i,:] = psi[K*(i+1)-1,:]

vecpsi = np.zeros(K*dis*(st+1))
for i in range(st+1): vecpsi[i*K*dis:(i+1)*K*dis] = psi[K*i:K*(i+1),:].T.ravel()

# Pi matrix (Hamilton 11.1.6)
piprime = np.zeros((K, K*p+dis*(q+1)))
for i in range(p): piprime[:, i*K:(i+1)*K] = phi[i*K:(i+1)*K,:]
for i in range(q+1): piprime[:, K*p+i*dis:K*p+(i+1)*dis] = theta[i*K:(i+1)*K,:]
pimat = piprime.T; smpi = pimat.ravel(order='F')

# Monte Carlo loop
mr_rnd = np.zeros(((st+1)*dis, monte)); ce = np.zeros((dis, monte))
for num in range(monte):
    smpi_rnd = np.random.multivariate_normal(smpi, np.kron(omega, iQ)/nobs)
    pi_rnd = smpi_rnd.reshape(K*p+dis*(q+1), K, order='F')
    psi_rnd = psimat(p,q,K,dis,st,pi_rnd)
    for i in range(dis):
        for j in range(st+1): mr_rnd[j+i*(st+1), num] = psi_rnd[K*(j+1)-1, i]
    for i in range(dis): ce[i,num] = mr_rnd[i*(st+1):i*(st+1)+4, num].sum()

# Quantiles
u95=round(monte*0.95)-1; l5=round(monte*0.05)-1
u90=round(monte*0.9)-1;  l10=round(monte*0.1)-1
mr_sort = np.sort(mr_rnd, axis=1); ce_sort = np.sort(ce, axis=1)
mr_low5  = mr_sort[:,l5].reshape(st+1,dis,order='F')
mr_low10 = mr_sort[:,l10].reshape(st+1,dis,order='F')
mr_up90  = mr_sort[:,u90].reshape(st+1,dis,order='F')
mr_up95  = mr_sort[:,u95].reshape(st+1,dis,order='F')
ce_sum = np.column_stack([ce_sort[:,l5], ce_sort[:,l10], ce_sort.mean(1), ce_sort[:,u90], ce_sort[:,u95]])

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

# Print table of mean responses (LSDV estimator)
ce_val = mr[0:4, :].sum(axis=0)

disaster_names = ['Droughts', 'Earthquakes', 'Floods', 'Storms']
disaster_cols  = [0, 1, 2, 3]

print('\n' + '='*72)
print(f'  Table: {dep} natural disaster shocks')
print('         (LSDV estimator)')
print('='*72)
print(f'{"Disaster":<16} {"Period":<20} {"Estimate":>14}')
print('-'*72)

for idx, j in enumerate(disaster_cols):
    name = disaster_names[idx]
    for i in range(4):
        sig = _sig_mr(mr[i,j], mr_up90[i,j], mr_low10[i,j], mr_up95[i,j], mr_low5[i,j])
        label = name if i == 0 else ''
        print(f'{label:<16} {"Year "+str(i):<20} {_fmt(mr[i,j], sig):>14}')
    sig_c = _sig_ce(ce_val[j], ce_sum[j,0], ce_sum[j,1], ce_sum[j,3], ce_sum[j,4])
    print(f'{"":<16} {"Cumulative effect":<20} {_fmt(ce_val[j], sig_c):>14}')

print('-'*72)
print(f'{"":>36} {"N countries":>14}  {N}')
print(f'{"":>36} {"Sample period":>14}  1960-2007')
print(f'{"":>36} {"Observations":>14}  {nobs}')
print('-'*72)
print('Note: * significant at one-tail 10% level, ** at 5% level.')
print('='*72)

elapsed = time.time() - t_start
print(f"\nInformation criteria (p={p}, q={q}):  AIC = {AIC:.4f}   SBC = {SBC:.4f}")
print(f"\nElapsed: {elapsed:.1f}s ({elapsed/60:.2f}min)")

# Graphs
x = np.arange(st+1); y = np.zeros(st+1)
for j, t in enumerate(['Drought','Earthquake','Flood','Storm']):
    plt.figure()
    plt.plot(x,mr[:,j],'-b',label='LSDV'); plt.plot(x,mr_up90[:,j],':b',label='90%')
    plt.plot(x,mr_low10[:,j],':b',label='10%'); plt.plot(x,y,'-k')
    plt.xlabel('Time (years)'); plt.title(f'{dep} {t} Shock'); plt.legend(); plt.tight_layout()
plt.show()
