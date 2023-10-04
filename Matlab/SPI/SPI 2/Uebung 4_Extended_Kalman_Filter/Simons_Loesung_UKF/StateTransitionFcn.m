function x_plus = StateTransitionFcn(x_minus, u, tSpan, p)
    odeFun = @(t,x) my_bioprocess_ode(t,x,u,p);
    [~,x] = ode45(odeFun,tSpan,x_minus);
    x_plus = x(end,:)';
end
