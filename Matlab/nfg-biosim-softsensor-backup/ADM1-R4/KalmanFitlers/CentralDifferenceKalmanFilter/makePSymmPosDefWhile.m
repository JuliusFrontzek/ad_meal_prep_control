function Pspd = makePSymmPosDefWhile(P,eps)
% slightly modifies a (possibly) non-symmetric, non-positive definite
% matrix P into a symmetric, possitive definite (spd) one (Pspd) that 
% enables cholesky decompsition
% - eps > 0: threshold of the smallest new eigenvalue

% 3 steps: 
% 1. make P symmetric: 
PTemp = (P + P')/2; 

% 2. increase negative eigenvalues, until smalles eigenvalue is above eps
myEps = eps;
counter = 0;
% also watch that eigenvalues remain real:
while (min(eig(PTemp)) < eps | isreal(eig(PTemp)) == false)
    
    [V,L] = eig(PTemp); % V - right eigenvectors, L - diag. with eigenvalues
    oldSmallestL = min(diag(L)) % print smallest eigenvalues
    lambdas = diag(L); 
    lambdas(lambdas < eps) = myEps; 

    % 3. create new (real) spd matrix through eigenvalue decomposition:
    V = real(V); 
    PTemp = V * real(diag(lambdas)) * V'; 
    
    % increase eps for next iteration
    myEps = 2*myEps

    % avoid infinite loop:
    counter = counter + 1
    if counter > 10
        break 
    end
end

Pspd = PTemp; 

newSmallesL = min(eig(Pspd)) % print smalles eigenvalue after eigen-decomposition

end