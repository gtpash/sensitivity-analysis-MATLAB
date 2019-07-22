% test SIR model Example 15.18 from [1]
%
% Author(s): Graham T. Pash | October 29, 2018
% 
% References
% [1] Smith, R. C. (2013). Uncertainty quantification: theory, implementation,
%     and applications (Vol. 12). SIAM.

clear; close all; clc

addpath('../')

% time vector for solution
tf = 5;
dt = 0.1;
t_vals = 0:dt:tf;

settings.model = @(x) SIRfun(t_vals,x);

% Large Degree of Interactions
alpha = 2; beta = 7; % k ~ Beta(alpha,beta)
settings.parameters(1) = struct('name', 'gamma', 'distpar', [.5 1],'dist','uniform');
settings.parameters(2) = struct('name', 'k', 'distpar',[alpha,beta],'dist','beta');
settings.parameters(3) = struct('name', 'r', 'distpar', [.5 1],'dist','uniform');
settings.parameters(4) = struct('name','delta','distpar', [.5 1],'dist','uniform');
settings.sample_points = 1e3;

output = saltelli_sobol(settings);

% Limited Interactions Case
alpha = 0.2; beta = 15;
settings.parameters(2) = struct('name', 'k', 'distpar',[alpha,beta],'dist','beta');
output2 = saltelli_sobol(settings);

%% auxillary functions
function dy = SIR_rhs(~,y,params)
S = y(1);
I = y(2);
R = y(3);
N = sum(y);

gamma = params(1);
k = params(2);
r = params(3);
d = params(4);

dy = [d*N - d*S - gamma*k*I*S;
    gamma*k*I*S - (r + d)*I;
    r*I - d*R];
end

function y = SIRfun(t,params)
y0 = [900;100;0];
[~,y] = ode45(@SIR_rhs,t,y0,'',params);
y = trapz(t, y(:,3)); % response is R(t,q) itegrated from 0 to 5
end