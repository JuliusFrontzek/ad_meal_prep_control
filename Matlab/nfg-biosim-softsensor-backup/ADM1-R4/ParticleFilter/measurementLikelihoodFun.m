function likelihood = measurementLikelihoodFun(particles,yMeas,pFix,R)
% computes the measurement likelihood function for each particle given the
% simulated outputs achieved with them
%
% Inputs:
%    particles - nParticles-by-nStates matrix that holds the particles
%
% Outputs:
%    likelihood - A vector with nParticles elements whose k-th
%                 element is the likelihood of the k-th particle

yPredicted = biogasmodell_mgl_ess_mat(particles,pFix);
[~,nMeasurements] = size(yPredicted); 

% Assume the RATIO of the measurement error to follow a Gaussian 
% distribution with zero mean, variance 0.2
% XY: hier noch die korrekten statistischen Daten eintragen:
mu = 0; % mean
% measVariance = R; % XY: hier evtl lieber direkt R einfügen!
invR = R \ eye(nMeasurements); 

% Use multivariate Gaussian probability density function, calculate
% likelihood of each particle
[nParticles,nStates] = size(particles);
likelihood = zeros(nParticles,1); % allocate memory
C = 1/sqrt(((2*pi)^nMeasurements*det(R)));   % Vorfaktor der mehr-dim. Gaußschen Normalverteilung 
% 
for kk = 1:nParticles
    % hier weicht Matlab von Dan Simon ab, der nicht durch yPredicted(kk)
    % teilt:
%     errorRatio = (yPredicted(kk) - yMeas)/yPredicted(kk);
    error = (yPredicted(kk,:) - yMeas); 
    v = error - mu;
    likelihood(kk) = C * exp(-0.5*(v*invR*v'));
end

end