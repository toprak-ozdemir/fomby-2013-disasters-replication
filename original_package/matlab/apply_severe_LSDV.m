%
%                        "apply_severe_LSDV.m                      
%  This is a program for deriving the LSDV estimates of the mean
%  responses of growth to natural disasters in the model of  
%
%       Thomas B. Fomby, Yuki Ikeda, and Norman V. Loayza, 
%       "The Growth Aftermath of Natural Disasters", Journal of Applied
%       Econometrics, forthcoming.  
% 
%  NOTE: 
%  (1) Severe disasters are differentiated.
%  (2) Select a file to load, in line 31.
%
%
%  Author: Yuki Ikeda
%  Email: yi22@georgetown.edu
%  Date: 8/10/2011
% ************************************************************************

clear all;
tic;

monte = 10000;                            % # of repetitions in the Monte Carlo loop
st = 6;                                   % # of periods for the forecast
dep = 'Mean Response of GDP growth to';   % GDP, agr., or non-agr. growth

% Set up parameter values      
p = 2;  q = 2;         % Panel VARX(p,q), 1<= p,q <= 3
   
% Read the data
vars = xlsread('dev_grow.xls');
[nobs,col] = size(vars);  % nobs = # of used observations

newid = vars(:,col);      % newid = country id (1~85)
N = max(newid);           % cross-sectional dimension
for i = 1:N;
    edge(i,1) = i;
end;
Ti = histc(newid,edge);  % Ti = time dimension of each country

acTi(1,1) = 0; 
for i = 1:N;             % acTi = accumulated Ti  
    acTi(i+1,1) = acTi(i) + Ti(i);
end;    

% Set model configuration
K = 6;       % # of endogenous vars
dis = 10;    % # of exogenous vars

Y = [vars(:,145) vars(:,94) vars(:,146:148) vars(:,2)]; 
LYs = [vars(:,149) vars(:,104) vars(:,150:152) vars(:,4) ...
       vars(:,153) vars(:,114) vars(:,154:156) vars(:,6)]; 
Xs = [vars(:,9:16) vars(:,1) vars(:,98) vars(:,17:24) vars(:,3) vars(:,108) ...
         vars(:,25:32) vars(:,5) vars(:,118)];  

     
% **********************************************************************************************         
nv = K*p + dis*(q+1);  % # of explanatory variables        
RH = [LYs Xs];         % Stacked right-hand-side vars 

% % Demean the series
dmY = zeros(nobs,K);
dmLYs = zeros(nobs,p*K);
dmXs = zeros(nobs,(p+1)*dis);

for nn = 1:N;
    dmY(1+acTi(nn):acTi(nn+1),:) = Y(1+acTi(nn):acTi(nn+1),:)-ones(Ti(nn),1)*mean(Y(1+acTi(nn):acTi(nn+1),:));
    dmLYs(1+acTi(nn):acTi(nn+1),:) = LYs(1+acTi(nn):acTi(nn+1),:)-ones(Ti(nn),1)*mean(LYs(1+acTi(nn):acTi(nn+1),:));
    dmXs(1+acTi(nn):acTi(nn+1),:) = Xs(1+acTi(nn):acTi(nn+1),:)-ones(Ti(nn),1)*mean(Xs(1+acTi(nn):acTi(nn+1),:));
end;

dmRH = [dmLYs dmXs]; % Demeaned tacked right-hand-side vars 

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

% Matrix inv(Q) 
iQ = inv((dmRH'*dmRH)/nobs);

% Estimated variance covariance matrix (omega hat)
eps = zeros(nobs,K);
eps = dmY - dmRH*deltahat;
omega = (eps'*eps)/nobs;
detomega = det((eps'*eps)/(nobs));

% Calculate matrix psi from deltahat
psi = psimat(p,q,K,dis,st,deltahat); 

% Set up phi matrix
phi = zeros(K*p,K);
for i = 1:p;    
    phi(1+K*(i-1):K*i,:) = deltahat(1+K*(i-1):K*i,:)';
end;

% Set up theta matrix
theta = zeros(K*(q+1),dis);
for i = 1:q+1;
    theta(1+(i-1)*K:K*i,:) = deltahat(1+K*p+dis*(i-1):K*p+dis*i,:)'; 
end;
    
%  Store mean responses
mr = zeros(st+1,dis);
for i = 1:st+1
        mr(i,:) = psi(K*i,:);  % Store only the effects on growth (the last one of end vars)
end;  

mr_sev = zeros(st+1,4);
for i = 1:4
        mr_sev(:,i) = mr(:,i) + mr(:,i+4);
end;  

% Vectorize psi matrix
vecpsi = zeros(K*dis*(st+1),1);
for i = 1:st+1;
    vecpsi(1+(i-1)*K*dis:K*dis*i) = reshape(psi(1+K*(i-1):K*i,:)',K*dis,1); 
end;    

% Set up Pi matrix (Hamilton 11.1.6)
piprime = zeros(K,K*p+dis*(q+1)); 
for i = 1:p;
    piprime(:,1+(i-1)*K:K*i) = phi(1+(i-1)*K:K*i,:);
end;
for i = 1:q+1;
    piprime(:,(K*p)+1+(i-1)*dis:(K*p)+i*dis) = theta(1+(i-1)*K:K*i,:);
end;    
pimat = piprime';
smpi = reshape(pimat,K*(K*p+dis*(q+1)),1);  % small pi (vector)

% **************************************************************
% Monte Carlo loop
mr_rnd = zeros((st+1)*dis,monte);
ce = zeros(dis,monte);  % Cumulative effects
mr_sev_rnd = zeros((st+1)*4,monte);
ce_sev = zeros(4,monte);  % Cumulative effects

for num = 1:monte;  % 'monte' times of iterations
    smpi_rnd = mvnrnd(smpi,kron(omega,iQ)/nobs);
    
    pi_rnd = reshape(smpi_rnd,K*p+dis*(q+1),K);    
    psi_rnd = psimat(p,q,K,dis,st,pi_rnd);

    % Store mean responses 
    for i = 1:dis
        for j = 1:st+1
            mr_rnd(j+(i-1)*(st+1),num) = psi_rnd(K*j,i);
        end;
    end;    
    
    % Store cumulative effects
    for i = 1:dis
        ce(i,num) = sum(mr_rnd(1+(i-1)*(st+1):4+(i-1)*(st+1),num));
    end;
   
    % Store mean responses 
    for i = 1:4;
        mr_sev_rnd(1+(i-1)*(st+1):i*(st+1),num) = mr_rnd(1+(i-1)*(st+1):i*(st+1),num) + mr_rnd(1+(4+i-1)*(st+1):(4+i)*(st+1),num);
    end;   

    % Store cumulative effects
    for i = 1:4;
        ce_sev(i,num) = sum(mr_sev_rnd(1+(i-1)*(st+1):4+(i-1)*(st+1),num));
    end; 
    
end;  % End of Monte Carlo loop
% **************************************************************

upper95 = round(monte*0.95); lower5 = round(monte*0.05);
upper90 = round(monte*0.9); lower10 = round(monte*0.1);

mr_sort = sort(mr_rnd,2);
ce_sort = sort(ce,2);

mr_low5 = reshape(mr_sort(:,lower5),st+1,dis);
mr_low10 = reshape(mr_sort(:,lower10),st+1,dis);
mr_up90 = reshape(mr_sort(:,upper90),st+1,dis);
mr_up95 = reshape(mr_sort(:,upper95),st+1,dis);

ce_sum = [ce_sort(:,lower5) ce_sort(:,lower10) mean(ce_sort,2) ce_sort(:,upper90) ce_sort(:,upper95)];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mr_sev_sort = sort(mr_sev_rnd,2);
ce_sev_sort = sort(ce_sev,2);

mr_sev_low5 = reshape(mr_sev_sort(:,lower5),st+1,4);
mr_sev_low10 = reshape(mr_sev_sort(:,lower10),st+1,4);
mr_sev_up90 = reshape(mr_sev_sort(:,upper90),st+1,4);
mr_sev_up95 = reshape(mr_sev_sort(:,upper95),st+1,4);

ce_sev_sum = [ce_sev_sort(:,lower5) ce_sev_sort(:,lower10) mean(ce_sev_sort,2) ce_sev_sort(:,upper90) ce_sev_sort(:,upper95)];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:4;
       disp('______');
    for i = 1:4;
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
    if ce_sum(j,2) > 0 & ce_sum(j,4) > 0 & ce_sum(j,1) < 0
        disp(' *');
    elseif ce_sum(j,1) > 0 & ce_sum(j,5) > 0
        disp('**');  
    elseif ce_sum(j,2) < 0 & ce_sum(j,4) < 0 & ce_sum(j,5) > 0
        disp(' *');
    elseif ce_sum(j,1) < 0 & ce_sum(j,5) < 0
        disp('**'); 
    elseif ce(j,3) > 0;
        disp('  ');
    else
        disp('  ');
    end;    
end;
disp('______');

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
    if ce_sev_sum(j,2) > 0 & ce_sev_sum(j,4) > 0 & ce_sev_sum(j,1) < 0
        disp(' *');
    elseif ce_sev_sum(j,1) > 0 & ce_sev_sum(j,5) > 0
        disp('**');    
    elseif ce_sev_sum(j,2) < 0 & ce_sev_sum(j,4) < 0 & ce_sev_sum(j,5) > 0
        disp(' *');
    elseif ce_sev_sum(j,1) < 0 & ce_sev_sum(j,5) < 0
        disp('**');    
    elseif ce_sev(j,3) > 0;
        disp('  ');
    else
        disp('  ');
    end;    
end;
disp('______');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(x,mr(:,5),'-b',x,mr_up90(:,5),':b',x,mr_low10(:,5),':b',x,y,'-k');
xlabel('Time (years)');
title([dep ' Severe Drought Shock ']);

figure;
plot(x,mr(:,6),'-b',x,mr_up90(:,6),':b',x,mr_low10(:,6),':b',x,y,'-k');
xlabel('Time (years)');
title([dep ' Severe Earthquake Shock ']);

figure; 
plot(x,mr(:,7),'-b',x,mr_up90(:,7),':b',x,mr_low10(:,7),':b',x,y,'-k');
xlabel('Time (years)');
title([dep ' Severe Flood Shock ']);

figure; 
plot(x,mr(:,8),'-b',x,mr_up90(:,8),':b',x,mr_low10(:,8),':b',x,y,'-k');
xlabel('Time (years)');
title([dep ' Severe Storm Shock ']);


toc;
disp(toc/60);


