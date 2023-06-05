%% Version
% (R2022b) Update 5
% Erstelldatum: 5.6.2023
% Autor: Simon Hellmann

function dxPdtNorm = dfP_dtNorm(xPNorm,uNorm,xiNorm,params,QNorm,fNorm,dfdxNorm,TxNum,TuNum)
                         % (xPNorm,u,xi,params,Q,f,dfdx)
% compute right-hand side of ODE of both states and state error covariance 
% matrix P. 
% dxPdt - right-hand side of ODE (of states and dynamics of state error
% covariance matrix)
% xP - states and covariance matrix (reshaped as vector), stacked above
% each other as column vector
% u - feedVolumeFlow
% xi - inlet concentrations (nStates,1)
% params - struct with stoichiometric coefficients a, aggregated constant
% parameters c and time-variant parameters th
% Q - power spectral density matrix of process noise
% f - function handle of ODEs of system equations
% dfdx - function handle of partial derivatives of df/dx 

    nStates = length(xiNorm); 
    
    % extract constant parameters out of struct: 
    th = params.th; 
    c = params.c; 
    a = params.a;
    
    dxPdtNorm = zeros(size(xPNorm));    % placeholder
    
    % separate states and covariance matrix:
    xNorm = xPNorm(1:nStates); 
    PNorm = reshape(xPNorm(nStates+1:end),[nStates,nStates]);
    
    % clipping of negative concentrations:
    % x(x < 0) = 0; 
    
    %% ODEs of states:
    dxdtNorm = fNorm(xNorm,uNorm,xiNorm,th,c,a,TxNum,TuNum);
    dxPdtNorm(1:nStates) = dxdtNorm; 
    
    %% ODEs of state error covariance matrix:
    
    % partial derivatives of the right-hand side w.r.t. states, evaluated
    % at current estimate x (which is actually xHat):
    FNorm = dfdxNorm(xNorm,uNorm,th,c,a,TxNum,TuNum);
%     Q = zeros(nStates);     % XY Rania
    dPdtNorm = FNorm*PNorm + PNorm*FNorm' + QNorm;  % dynamics of state error covariance matrix
    % reshape matrix as long column vector and append values dxPdt:
    dxPdtNorm(nStates+1:end) = reshape(dPdtNorm,[],1);

end