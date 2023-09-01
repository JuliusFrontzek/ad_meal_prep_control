function jacobian = dfdp(x,a)
% evaluate the jacobian of f w.r.t. the parameters p at current state x
% a - array of stoichiometric constants acc. to petersen matrix

jacobian = [a(1,1).*x(3),a(1,2).*x(4),a(1,3).*x(5),0;
            a(2,1).*x(3),a(2,2).*x(4),a(2,3).*x(5),0;
            -x(3),0,0,a(5,4).*x(6);
            0,-x(4),0,a(6,4).*x(6);
            0,0,-x(5),a(7,4).*x(6);
            a(8,1).*x(3),a(8,2).*x(4),a(8,3).*x(5),-x(6);
            0,0,0,0;
            0,0,0,0];

end