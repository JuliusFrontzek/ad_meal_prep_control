function dfdp = dfdp(x,u,p,c,a)
% dfdp - partial derivative of ode w.r.t. parameters, evaluated at current
% state x
% x - state variables
% u - feed volume flow
% a - stoichiometric constants from Petersen Matrix
% c - summarized model parameters

dfdp = [a(1,1).*x(5),a(2,1).*x(6),a(1,3).*x(7),0;
        a(2,1).*x(5),a(2,2).*x(6),a(2,3).*x(7),0;
        -a(3,1).*x(5),a(3,2).*x(6),-a(3,3).*x(7),0;
        -a(4,1).*x(5),-a(4,2).*x(6),-a(4,3).*x(7),0;
        -x(5),0,0,a(5,4).*x(8);
        0,-x(6),0,a(6,4).*x(8);
        0,0,-x(7),a(7,4).*x(8);
        a(8,1).*x(5),a(8,2).*x(6),a(8,3).*x(7),-x(8);
        0,0,0,0;
        0,0,0,0;
        0,0,0,0];

end