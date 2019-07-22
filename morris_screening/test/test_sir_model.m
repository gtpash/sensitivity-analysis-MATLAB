% test morris screening
% Example 15.18 from [1], SIR Model
% 
% Author(s): Graham T. Pash | October 23, 2018

% References
% [1] Smith, R. C. (2013). Uncertainty quantification: theory, implementation,
%     and applications (Vol. 12). SIAM.

clear; close all; clc

addpath('../')

alpha = .2; beta = 15; % k ~ Beta(alpha,beta)

% set up parameters
settings.parameters(1) = struct('name', 'gamma', 'dist','uniform','distpar',[0.5 1]);
settings.parameters(2) = struct('name', 'k', 'dist','beta','distpar',[alpha,beta]);
settings.parameters(3) = struct('name', 'r', 'dist','uniform','distpar',[0.5 1]);
settings.parameters(4) = struct('name','delta','dist','uniform','distpar',[0.5 1]);

% model settings
settings.sample_points = 1e3;
settings.use_method = 'fd';
settings.delta = 1e-6;

tf = 5;
dt = 0.1;
t_vals = 0:dt:tf;
settings.model = @(x) SIRfun(t_vals,x);

out = morris_screening(settings);

%% auxillary functions
function dy = SIRsys(~,y,params)
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
[~,y] = ode45(@SIRsys,t,y0,'',params);
y = trapz(t, y(:,3));
end