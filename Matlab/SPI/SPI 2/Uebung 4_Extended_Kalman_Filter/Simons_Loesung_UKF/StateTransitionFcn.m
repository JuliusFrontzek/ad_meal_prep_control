function xPlus = StateTransitionFcn(xMinus, u, tSpan, p)
    odeFun = @(t,x) my_bioprocess_ode(t,x,u,p);
    [~,x] = ode45(odeFun,tSpan,xMinus);
    xPlus = x(end,:)';
end
