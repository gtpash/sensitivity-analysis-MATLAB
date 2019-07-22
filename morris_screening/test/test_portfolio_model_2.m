% test portfolio model with sigma-scaled elementary effect
% Example 15.16 from [1]
% 
% Author(s): Graham T. Pash | October 23, 2018

% References
% [1] Smith, R. C. (2013). Uncertainty quantification: theory, implementation,
%     and applications (Vol. 12). SIAM.

clear; close all; clc;

addpath('../')

c1 = 2;
c2 = 1;

% obtain standard deviations through MCMC analysis, such as DRAM
sig1 = 1;
sig2 = 3;
sigy = c1^2*sig1^2 + c2^2*sig2^2; % model variance

settings.model = @(x) c1*x(1) + c2*x(2);

settings.parameters(1) = struct('name', 'Q1', 'distpar', [0, sig1^2], 'dist', 'normal');
settings.parameters(2) = struct('name', 'Q2', 'distpar', [0, sig2^2], 'dist', 'normal');
settings.sample_points = 200;
settings.delta = 1e-6;
settings.use_method = 'complex';
settings.observation_error = sigy;

output = morris_screening(settings);

fprintf('mu_i_star = [%4.4f, %4.4f]\nsigma_i = [%4.2e, %4.2e]\n', output.mu_i_star, output.sigma_i);
