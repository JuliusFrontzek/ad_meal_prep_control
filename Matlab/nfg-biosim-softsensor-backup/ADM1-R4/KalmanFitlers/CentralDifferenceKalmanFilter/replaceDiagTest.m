% replace the diagonal of one Matrix A through a vector diagSoll: 

A = magic(3);   
diagSoll = ones(3,1); 

ANoDiag = A - diag(diag(A)); 
ANewDiag = ANoDiag + diag(diagSoll); 


%% does any arbitrary matrix become symmetric if it is multiplied with the 
% same matrix (untransposed & transposed) from both sides?

B = magic(3); 
C = magic(3); 

D = C*B*C'; 



