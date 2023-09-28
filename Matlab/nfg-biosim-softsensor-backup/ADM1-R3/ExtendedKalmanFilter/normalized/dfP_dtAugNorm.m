%% Version
% (R2022b) Update 5
% Erstelldatum: 01.08.2023
% Autor: Simon Hellmann

function dxPNormdt = dfP_dtAugNorm(xPNorm,uNorm,xiNorm,params,QNorm, ...
                                   fNorm,dfdxNorm,TxNum,TuNum,...
                                   nStates,flagAugmented)
% compute right-hand side of ODEs of both normalized (norm.) states and 
% norm. state error covariance matrix P. Respect both cases (augmented and
% non-augmented)

% dxPNormdt - right-hand side of ODE (of norm. states and dynamics of norm. 
% state error covariance matrix)
% xPNorm - norm. states and covariance matrix (reshaped as vector), stacked above
% each other as column vector
% uNorm - norm. feedVolumeFlow
% xiNorm - norm. inlet concentrations (nStates,1)
% params - struct with stoichiometric coefficients a, aggregated constant
% parameters c and time-variant parameters th
% QNorm - power spectral density matrix of norm. process noise
% fNorm - function handle of norm. ODEs of system equations
% dfdxNorm - function handle of norm. partial derivatives of df/dx 
% TxNum, TuNum - normalization vectors of states and inputs
% nStates - # states (without sample state augmentation)
% flagAugmented - 0: non-augmented, 1: augmented

    % extract constant parameters out of struct: 
    th = params.th; 
    c = params.c; 
    a = params.a;
    
    dxPNormdt = zeros(size(xPNorm));    % allocate memory 
    
    % separate states and covariance matrix:
    if flagAugmented == 1
        xAugNorm = xPNorm(1:2*nStates); 
        PAugNorm = reshape(xPNorm(2*nStates+1:end),[2*nStates,2*nStates]);
        xNorm = xAugNorm(1:nStates); % ignore the sample state for further computations
        PNorm = PAugNorm; 
    else
        xNorm = xPNorm(1:nStates); 
        PNorm = reshape(xPNorm(nStates+1:end),[nStates,nStates]);
    end    
    
    % clipping of negative concentrations:
    % x(x < 0) = 0; 
    
    %% ODEs of states:
    dxdtNorm = fNorm(xNorm,uNorm,xiNorm,th,c,a,TxNum,TuNum);                

    if flagAugmented == 1
        dxsNormdt = zeros(nStates,1); % maintain sample state
        dxAugNormdt = [dxdtNorm;dxsNormdt]; 
        dxPNormdt(1:2*nStates) = dxAugNormdt; 
    else 
        dxPNormdt(1:nStates) = dxdtNorm; 
    end
    
    %% ODEs of state error covariance matrix:
    
    % partial derivatives of the right-hand side w.r.t. states, evaluated
    % at current estimate x (which is actually xHat):
    FNorm = dfdxNorm(xNorm,uNorm,th,c,a,TxNum,TuNum);            
    if flagAugmented == 1
        % augment F- and Q-matrix:
        FNorms = zeros(nStates);    % for normalized sample-states
        FAugNorm = blkdiag(FNorm,FNorms);  % augmented, normalized F-matrix
%         Q = zeros(nStates);   % XY Rania
        QNorms = zeros(nStates);    % for sample-states
        QAugNorm = blkdiag(QNorm,QNorms);     % augmented Q-matrix
        
        dPAugNormdt = FAugNorm*PAugNorm + PAugNorm*FAugNorm' + QAugNorm;  % dynamics of normalized state error covariance matrix
        % reshape matrix as long column vector and append values dxPdt:
        dxPNormdt(2*nStates+1:end) = reshape(dPAugNormdt,[],1);
    else 
%         Q = zeros(nStates);   % XY Rania
        dPNormdt = FNorm*PNorm + PNorm*FNorm' + QNorm;  % dynamics of state error covariance matrix
        % reshape matrix as long column vector and append values dxPdt:
        dxPNormdt(nStates+1:end) = reshape(dPNormdt,[],1);
    end
    
end