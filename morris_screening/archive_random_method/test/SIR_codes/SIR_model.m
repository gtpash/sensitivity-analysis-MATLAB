function y = SIR_model(N, t_vals, Y0, ode_options, theta)

params.theta.gamma = theta(1);
params.theta.k = theta(2);
params.theta.r = theta(3);
params.theta.delta = theta(4);
params.N = N;

[t,Y] = ode45(@SIR_rhs,t_vals,Y0, ode_options, params);
Y(:,3) = N - Y(:,1) - Y(:,2);

y = trapz(t, Y(:,3));