function dhdx = dhdx(x)

%% Zustände
mX = x(1);
mS = x(2);
V = x(3);

%% Ableitung
  
dhdx = [1/V 0 -mX/V^2; 
        0   0 1];
end
