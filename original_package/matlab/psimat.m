% 
%                        "psimat.m" (function)        
%  A function for calculating matrix psi from phi and theta in
%  apply_BSBC.m, apply_severe_BSBC.m, and meanpsi.m.
%         
%
%  Author: Yuki Ikeda
%  Email: yi22@georgetown.edu
%  Date: 8/10/2011
% ************************************************************************

function psi = psimat(p,q,nen,dis,st,beta)


% Set up phi matrix    
phi = zeros(nen*p,nen);
for i = 1:p;    
    phi(1+nen*(i-1):nen*i,:) = beta(1+nen*(i-1):nen*i,:)';
end;
phi1 = phi(1:nen,:); 
phi2 = phi(1+nen:2*nen,:); phi3 = zeros(nen,nen);

% Set up theta matrix
theta = zeros(nen*(q+1),dis);
for i = 1:q+1;
    theta(1+(i-1)*nen:nen*i,:) = beta(1+nen*p+dis*(i-1):nen*p+dis*i,:)'; 
end;
theta0 = theta(1:nen,:); theta1 = theta(1+nen:2*nen,:); 
theta2 = theta(1+2*nen:3*nen,:); theta3 = zeros(nen,dis);

% Calculate psi matrix
psi = zeros(nen*(st+1),dis); 
psi(1:nen,:) = theta0;  % psi0 = theta0 
psi(1+nen:2*nen,:) = theta1 + phi1*psi(1:nen,:);  % psi1 =  theta1 + phi1*psi0
psi(1+2*nen:3*nen,:) = theta2 + phi1*psi(1+nen:2*nen,:) + phi2*psi(1:nen,:);  % psi2 = theta2 + phi1*psi1 + phi2*psi0
psi(1+3*nen:4*nen,:) = theta3 + phi1*psi(1+2*nen:3*nen,:) + phi2*psi(1+nen:2*nen,:) + phi3*psi(1:nen,:);  % psi3

for i = 4:st;
    psi(1+nen*i:nen*(i+1),:) = phi1*psi(1+nen*(i-1):nen*i,:)+ ... 
                           phi2*psi(1+nen*(i-2):nen*(i-1),:) + ...
                           phi3*psi(1+nen*(i-3):nen*(i-2),:);  
end;
