function Q = updateQ(x,u,AC,kQ)
% compute time-varianz spectral density matrix Q acc. to Schneider &
% Georgakis (2013)

    Cov_pp = AC.Cov_pp;
    p = AC.th; 
    a = AC.a; 
    c = AC.c;
    Jp = dfdp(x,u,p,c,a);   % jacobian of ode w.r.t. identified paramaters p
    Q = kQ*Jp*Cov_pp*Jp'; 

end