function dxdt = my_bioprocess_ode(t,x,u,p)
    
    dxdt = nan(3,1);    % allocate memory
    
    % Parameter:
    mumax = p(1);
    KS = p(2);
    YXS = p(3);
       
    % Grenzen der Zustände:hard reset nedativer Konzentrationen:
    for idxX = 1:length(x)
	    if x(idxX) < 0
		    x(idxX) = 0;
	    end % if
    end % for idxX
    
    % Definition der Zustände:
    mX = x(1);
    mS = x(2);
    V = x(3);
    
    % Definition der Stellgrößen:
    q = u;
    cSF = 100;
    
    % Zwischengrößen:
    cS = mS/V;
    mu = mumax * cS / (cS + KS);
    
    % DGLs für Zustände:
    dxdt(1) = mu * mX;
    dxdt(2) = -1/YXS * mu * mX + cSF * q;
    dxdt(3) = q;

end