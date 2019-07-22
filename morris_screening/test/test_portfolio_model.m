% test morris screening
%
% Author(s): Paul R. Miles | August 8, 2018

clear; close all; clc;

addpath('../')

c1 = 2;
c2 = 1;

settings.model = @(x) c1*x(1) + c2*x(2);

settings.parameters(1) = struct('name', 'Q1', 'distpar', [0, 1], 'dist', 'normal');
settings.parameters(2) = struct('name', 'Q2', 'distpar', [0, 9], 'dist', 'normal');
settings.sample_points = 200;
settings.delta = 1e-6;
settings.use_method = 'complex';

output = morris_screening(settings);

fprintf('mu_i_star = [%4.4f, %4.4f]\nsigma_i = [%4.2e, %4.2e]\n', output.mu_i_star, output.sigma_i);
