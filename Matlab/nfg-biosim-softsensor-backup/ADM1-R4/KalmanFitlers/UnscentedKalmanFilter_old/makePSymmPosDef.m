function Pspd = makePSymmPosDef(P,eps)
% slightly modifies a (possibly) non-symmetric, non-positive definite
% matrix P into a symmetric, possitive definite (spd) one (Pspd) that 
% enables cholesky decompsition
% - eps > 0: threshold of the smallest new eigenvalue

% 3 steps: 
% 1. make P symmetric: 
PTemp = (P + P')/2; 

% 2. increase negative eigenvalues up to eps
[V,L] = eig(PTemp); % V - right eigenvectors, L - diag. with eigenvalues
% L = diag(lam1, lam2, ..., lamn)
L(L < eps) = eps; 

% 3. create new spd matrix through eigenvalue decomposition:
Pspd = V * L * V'; 

end