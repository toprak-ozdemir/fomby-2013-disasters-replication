% 
%                        "btst.m" (function)        
%  A function for calculating deltatil_b in apply_BSBC.m and
%  apply_severe_BSBC.m.
%                                              
%
%  Author: Yuki Ikeda
%  Email: yi22@georgetown.edu
%  Date: 8/10/2011
% ************************************************************************

function deltatil_b = btst(B,nobs,nv,p,K,N,A,Y,LYs,Xs,Ti,acTi,D,deltahat)


% Set initial guess of deltatil = original LSDV estimate 
RH = [LYs Xs];
deltatil_b = zeros(nv,K*B);     % storage for delta tilda b (bootstrap samples)
deltatil_old = deltahat;        % initial guess of delta tilda

% Given deltatil, calculate etatil, epstil, and epstil_r 

alphatil_pre = D'*(Y-RH*deltatil_old);
alphatil = zeros(N,K);
for i = 1:N;
    alphatil(i,:) = alphatil_pre(i,:)/Ti(i);
end;    
D_alphatil = D*alphatil;

epstil = Y - RH*deltatil_old - D*alphatil;

% Rescaling residuals
M = diag((eye(nobs) - A*RH*inv(RH'*A*RH)*RH')*A);   
eom = zeros(nobs,K);
for i = 1:K;
    eom(:,i) = epstil(:,i)./(M.^0.5);    % eps over sqrt(diag(M))
end;

epstil_r = zeros(nobs,K);
for i = 1:N;
    epstil_r(1+acTi(i):acTi(i+1),:) = ((Ti(i)/(Ti(i)-1))^0.5)*...
        (eom(1+acTi(i):acTi(i+1),:) - ones(Ti(i),1)*mean(eom(1+acTi(i):acTi(i+1),:)));
end;

% Bootstrap sampling (B times)
for brep = 1:B; 

    % Resampling
    rb = ceil(rand(nobs,1)*nobs);
    for i = 1:nobs;
        epstil_b(i,:) = epstil_r(rb(i),:);
    end;    

    % Set ititial lagged values
    ly_b = zeros(nobs,p*K);
    for i = 1:N;
        ly_b(acTi(i)+1,:) = LYs(acTi(i)+1,:);
        ly_b(acTi(i)+2,1+K:2*K) = LYs(acTi(i)+2,1+K:2*K);
    end;        

    % Set y 
    y_b = zeros(nobs,K);

    for nn = 1:N;
        for i = 1+acTi(nn):acTi(nn+1);
            % if 1 <= t <= Ti-2
            if (i >= 1+acTi(nn)) & (i <= acTi(nn+1)-2);
                y_b(i,:) = [ly_b(i,:) Xs(i,:)]*deltatil_old ...
                     + D_alphatil(i,:) + epstil_b(i,:);
                ly_b(i+1,1:K) = y_b(i,:);  % ly(t+1) = y(t)
                ly_b(i+2,1+K:2*K) = y_b(i,:);  % l2y(t+2) = y(t)

            % if t = Ti-1     
            elseif i == acTi(nn+1)-1;    
                y_b(i,:) = [ly_b(i,:) Xs(i,:)]*deltatil_old ...
                    + D_alphatil(i,:) + epstil_b(i,:);
                ly_b(i+1,1:K) = y_b(i,:);  % ly(t+1) = y(t)             

            % if t = Ti
            elseif i == acTi(nn+1);
                y_b(i,:) = [ly_b(i,:) Xs(i,:)]*deltatil_old ...
                    + D_alphatil(i,:) + epstil_b(i,:);
            else;
            end;    
        end;    
    end;       

    % LSDV estimator   
    W_b = [ly_b Xs];
    deltatil_b(:,1+K*(brep-1):K*brep) = inv(W_b'*A*W_b)*W_b'*A*y_b;

end;    % end of brep (bootstrap)

