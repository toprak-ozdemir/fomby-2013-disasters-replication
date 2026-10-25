% 
%                        "meanpsi.m" (funtion)      
%  A function for taking the mean of psi_b over B times replications in
%  apply_BSBC.m and apply_severe_BSBC.m.
%                                        
%
%  Author: Yuki Ikeda
%  Email: yi22@georgetown.edu
%  Date: 8/10/2011
% ************************************************************************

function m_psi_b = meanpsi(p,q,nen,dis,st,K,B,deltatil_b)

m_psi_b = zeros(nen*(st+1),dis); 

for i = 1:B;
    m_psi_b = m_psi_b + psimat(p,q,nen,dis,st,deltatil_b(:,(i-1)*K+1:K*i));
end;

m_psi_b = m_psi_b/B;