% test landau energy example from Spring 2016 MA540
%
% Author(s): Graham T. Pash | October 23, 2018

clear; close all; clc;

addpath('../')

P = linspace(0,.8);

% analytically integrate to obtain scalar response
settings.model = @(x) sum(.8.^[3 5 7]./[3 5 7].*x');

settings.parameters(1) = struct('name', 'a1', 'distpar', [-389.4, 0.2], 'dist', 'uniform');
settings.parameters(2) = struct('name', 'a11', 'distpar', [761.3, 10], 'dist', 'normal');
settings.parameters(3) = struct('name', 'a111', 'distpar', [61.46, 0.2], 'dist', 'uniform');
settings.sample_points = 2000;
settings.delta = 1e-6;
settings.use_method = 'cs';
settings.scale_grad = true;

output = morris_screening(settings);

fprintf('mu_i_star = [%4.4f, %4.4f, %4.4f]\nsigma_i = [%4.2e, %4.2e, %4.2e]\n', output.mu_i_star, output.sigma_i);
