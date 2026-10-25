% 
%                        "apply_severe_BSBC.m"           
%  This is a program for deriving the BSBC estimates of the mean
%  responses of growth to natural disasters in the model of
%
%       Thomas B. Fomby, Yuki Ikeda, and Norman V. Loayza, 
%       "The Growth Aftermath of Natural Disasters", Journal of Applied
%       Econometrics, forthcoming.  
% 
%  NOTE:
%  (1) Severe disasters are differentiated.
%  (2) Select a file to load, in line 38.
%  (3) Set model configuration, in line 54-63.
%     Y: Endogenous variables (no lag)
%     LYs: [Endo.variables [lag1], Endo. variables [lag2], Endo.variables[lag3]]
%     Xs: [Exo.variables [lag0], Exo.variables [lag1], Exo.variables[lag2], Exo.variables[lag3]]
%
%     For LYs and Xs, the order is as: [varA[lag t], varB[lag t], varC[lag t], varA[lag t+1], varB[lag t+1], ...]
%     Not  [varA[lag t], varA[lag t+1], varA[lag t+2], varB[lag t], ...]
%
%
%  Author: Yuki Ikeda
%  Email: yi22@georgetown.edu
%  Date: 8/10/2011
% *************************************************************************************************************

clear all;
tic;

% Set parameter values
B = 1000;              % # of bootstrap samples     
st = 6;                % # of periods for the forecast
dep = 'Mean Response of GDP growth to';   % GDP, agr.growth, or nonagr.growth

p = 2;  q = 2;         % Panel VARX(p,q), 1<= p,q <= 3

% Load data
vars = xlsread('dev_grow.xls');

[nobs,col] = size(vars); % nobs = # of used observations

newid = vars(:,col);     % newid = country id 
N = max(newid);          % cross-sectional dimension
for i = 1:N;
    edge(i,1) = i;
end;
Ti = histc(newid,edge);  % Ti = time dimension of each country

acTi(1,1) = 0; 
for i = 1:N;             % acTi = accumulated Ti  
    acTi(i+1,1) = acTi(i) + Ti(i);
end;    

% Set model configuration
K = 6;                   % # of endogenous vars
dis = 10;                % # of exogenous vars (6 + 4 severe)

% % Benchmark model [Table V in the paper]
Y = [vars(:,145) vars(:,94) vars(:,146:148) vars(:,2)]; 
LYs = [vars(:,149) vars(:,104) vars(:,150:152) vars(:,4) ...
       vars(:,153) vars(:,114) vars(:,154:156) vars(:,6)];
Xs = [vars(:,9:16) vars(:,1) vars(:,98) vars(:,17:24) vars(:,3) vars(:,108) ...
      vars(:,25:32) vars(:,5) vars(:,118)];   
  
% **********************************************************************************************************      
nv = K*p + dis*(q+1);    % # of explanatory variables per equation
RH = [LYs Xs];           % Stacked right-hand-side vars

% % Demean the series
% dmY = zeros(nobs,K);
% dmLYs = zeros(nobs,p*K);
% dmXs = zeros(nobs,(p+1)*dis);
% 
% for nn = 1:N;
%     dmY(1+acTi(nn):acTi(nn+1),:) = Y(1+acTi(nn):acTi(nn+1),:)-ones(Ti(nn),1)*mean(Y(1+acTi(nn):acTi(nn+1),:));
%     dmLYs(1+acTi(nn):acTi(nn+1),:) = LYs(1+acTi(nn):acTi(nn+1),:)-ones(Ti(nn),1)*mean(LYs(1+acTi(nn):acTi(nn+1),:));
%     dmXs(1+acTi(nn):acTi(nn+1),:) = Xs(1+acTi(nn):acTi(nn+1),:)-ones(Ti(nn),1)*mean(Xs(1+acTi(nn):acTi(nn+1),:));
% end;
% 
% dmRH = [dmLYs dmXs];


% Derive LSDV estimate of original series

% Matrix D
D = zeros(nobs,N);
for nn = 1:N;
    for i = 1+acTi(nn):acTi(nn+1);
            D(i,nn) = 1;
    end;
end;    

% Matrix A 
A = zeros(nobs,nobs);
for i = 1:N;
    A(acTi(i)+1:acTi(i+1),acTi(i)+1:acTi(i+1)) = ...
        eye(Ti(i)) - (1/Ti(i))*ones(Ti(i),1)*ones(1,Ti(i)); 
end;  

% LSDV estimator
deltahat = inv(RH'*A*RH)*RH'*A*Y;

% To check deltahat = demeaned OLS
% deltahat_dmtest = inv(dmRH'*dmRH)*dmRH'*dmY;  


% ********************* %
%  Bootstrap procedure
% ********************* %  

% Derive the mean of psi_b, using subroutine "btst"
deltatil_b = btst(B,nobs,nv,p,K,N,A,Y,LYs,Xs,Ti,acTi,D,deltahat);
m_psi_b = meanpsi(p,q,K,dis,st,K,B,deltatil_b);

% Calculate psi of Pesaran and Zhao 
psi_pz = 2*psimat(p,q,K,dis,st,deltahat) - m_psi_b;

% Store mean responses
mr = zeros(st+1,8);
for i = 1:st+1;
        mr(i,:) = psi_pz(K*i,1:8);  % Store only the effects on growth (the last one of end vars)
end;  

mr_sev = zeros(st+1,4);
for i = 1:4
        mr_sev(:,i) = mr(:,i) + mr(:,i+4);
end;  

% Derive BSBC estimate for each of deltatil_b
mr_store = zeros(st+1,8*B);
mr_sev_store = zeros(st+1,4*B);

for h = 1:B;
    beta_h = deltatil_b(:,(h-1)*K+1:K*h);
    deltatil_b_h = btst(B,nobs,nv,p,K,N,A,Y,LYs,Xs,Ti,acTi,D,beta_h);
    m_psi_b_h =  meanpsi(p,q,K,dis,st,K,B,deltatil_b_h);
    psi_pz_h = 2*psimat(p,q,K,dis,st,beta_h) - m_psi_b_h;
    
    for i = 1:st+1;
        mr_store(i,(h-1)*8+1:8*h) = psi_pz_h(K*i,1:8);  % Store only the effects on growth (the last one of end vars)
        mr_sev_store(i,(h-1)*4+1:4*h) = mr_store(i,(h-1)*8+1:(h-1)*8+4) + mr_store(i,(h-1)*8+5:h*8);
    end;     
    disp(h);    
end;   

mr_store_long = zeros(B,8*(st+1));
mr_sev_store_long = zeros(B,4*(st+1));

for h = 1:B;
	mr_store_long(h,:) = reshape(mr_store(:,1+8*(h-1):8*h)',1,8*(st+1));
	mr_sev_store_long(h,:) = reshape(mr_sev_store(:,1+4*(h-1):4*h)',1,4*(st+1));
end;

% Sort mean responses
mr_store_long = sort(mr_store_long,1);
mr_low5 = reshape(mr_store_long(round(B*0.05),:),8,st+1)';
mr_low10 = reshape(mr_store_long(round(B*0.1),:),8,st+1)';
mr_up90 = reshape(mr_store_long(round(B*0.90),:),8,st+1)';
mr_up95 = reshape(mr_store_long(round(B*0.95),:),8,st+1)';

mr_sev_store_long = sort(mr_sev_store_long,1);
mr_sev_low5 = reshape(mr_sev_store_long(round(B*0.05),:),4,st+1)';
mr_sev_low10 = reshape(mr_sev_store_long(round(B*0.1),:),4,st+1)';
mr_sev_up90 = reshape(mr_sev_store_long(round(B*0.90),:),4,st+1)';
mr_sev_up95 = reshape(mr_sev_store_long(round(B*0.95),:),4,st+1)';


% Store cumulative effects
ce_store = sum(mr_store(1:4,:));
ce_store_long = zeros(B,8);

ce_sev_store = sum(mr_sev_store(1:4,:));
ce_sev_store_long = zeros(B,4);

for h = 1:B;
	ce_store_long(h,:) = ce_store(:,1+8*(h-1):8*h);
   	ce_sev_store_long(h,:) = ce_sev_store(:,1+4*(h-1):4*h);
end;
ce_store_long = sort(ce_store_long,1);
ce_sev_store_long = sort(ce_sev_store_long,1);

ce = zeros(5,8);
ce(1,:) = ce_store_long(round(B*0.05),:);
ce(2,:) = ce_store_long(round(B*0.1),:);
ce(3,:) = sum(mr(1:4,:));
ce(4,:) = ce_store_long(round(B*0.9),:);
ce(5,:) = ce_store_long(round(B*0.95),:);

ce_sev = zeros(5,4);
ce_sev(1,:) = ce_sev_store_long(round(B*0.05),:);
ce_sev(2,:) = ce_sev_store_long(round(B*0.1),:);
ce_sev(3,:) = sum(mr_sev(1:4,:));
ce_sev(4,:) = ce_sev_store_long(round(B*0.9),:);
ce_sev(5,:) = ce_sev_store_long(round(B*0.95),:);


for j = 1:4;      % 4 disasters
       disp('______');

    for i = 1:4;  % 0~3 periods
        if mr_up90(i,j) > 0 & mr_low10(i,j) > 0 & mr_low5(i,j) < 0
            disp(' *');
        elseif mr_up95(i,j) > 0 & mr_low5(i,j) > 0 
            disp('**');
        elseif  mr_up90(i,j) < 0 & mr_low10(i,j) < 0 & mr_up95(i,j) > 0
            disp(' *');
        elseif mr_up95(i,j) < 0 & mr_low5(i,j) < 0
            disp('**');    
        elseif mr(i,j) > 0
            disp('  ');
        else
            disp('  ');
        end;    
    end;    
    if ce(4,j) > 0 & ce(2,j) > 0 & ce(1,j) < 0
        disp(' *');
    elseif ce(5,j) > 0 & ce(1,j) > 0
        disp('**');  
    elseif ce(4,j) < 0 & ce(2,j) < 0 & ce(5,j) > 0
        disp(' *');
    elseif ce(5,j) < 0 & ce(1,j) < 0
        disp('**'); 
    elseif ce(3,j) > 0;
        disp('  ');
    else
        disp('  ');
    end;    

end;
disp('______');

for j = 1:4;
    disp('______');
    for i = 1:4;
        if mr_sev_up90(i,j) > 0 & mr_sev_low10(i,j) > 0 & mr_sev_low5(i,j) < 0
            disp(' *');
        elseif mr_sev_up95(i,j) > 0 & mr_sev_low5(i,j) > 0 
            disp('**')
        elseif  mr_sev_up90(i,j) < 0 & mr_sev_low10(i,j) < 0 & mr_sev_up95(i,j) > 0
            disp(' *');
        elseif mr_sev_up95(i,j) < 0 & mr_sev_low5(i,j) < 0
            disp('**');        
        elseif mr_sev(i,j) > 0
            disp('  ');
        else
            disp('  ');
        end;    
    end;    
    if ce_sev(4,j) > 0 & ce_sev(2,j) > 0 & ce_sev(1,j) < 0
        disp(' *');
    elseif ce_sev(5,j) > 0 & ce_sev(1,j) > 0
        disp('**');  
    elseif ce_sev(4,j) < 0 & ce_sev(2,j) < 0 & ce_sev(5,j) > 0
        disp(' *');
    elseif ce_sev(5,j) < 0 & ce_sev(1,j) < 0
        disp('**'); 
    elseif ce_sev(3,j) > 0;
        disp('  ');
    else
        disp('  ');
    end;   
end;
disp('______');


save dev_grow_sev_1000.mat
toc;
disp(toc/60);


% Draw the graphs   
x = (0:1:st)'; y = zeros(st+1,1);
figure;
plot(x,mr(:,1),'-b',x,mr_up90(:,1),':b',x,mr_low10(:,1),':b',x,y,'-k');
xlabel('Time (years)');
title([dep ' Drought Shock ']);

figure;
plot(x,mr(:,2),'-b',x,mr_up90(:,2),':b',x,mr_low10(:,2),':b',x,y,'-k');
xlabel('Time (years)');
title([dep ' Earthquake Shock ']);

figure; 
plot(x,mr(:,3),'-b',x,mr_up90(:,3),':b',x,mr_low10(:,3),':b',x,y,'-k');
xlabel('Time (years)');
title([dep ' Flood Shock ']);

figure; 
plot(x,mr(:,4),'-b',x,mr_up90(:,4),':b',x,mr_low10(:,4),':b',x,y,'-k');
xlabel('Time (years)');
title([dep ' Storm Shock ']);
