function f = odeFunSymTrial(x, u, a)

    % dynamic equations:
    f = [-a(2,1)*x(1)*u - a(1,1)*x(2);        
         -a(1,2)*x(3) - a(2,2)*u*x(2)^2; 
         -a(3,1)*x(2) - u^2*a(3,2)];

end 
