% test landau energy model from MA540 
%
% Author(s): Graham T. Pash | October 23, 2018

clear; close all; clc;

addpath('../')

P = linspace(0,.8);
% analytically integrate to obtain scalar response
settings.model = @(x) sum(.8.^[3 5 7]./[3 5 7].*x);

settings.parameters(1) = struct('name', 'a1', 'distpar', [-389.4, 0.2], 'dist', 'uniform');
settings.parameters(2) = struct('name', 'a11', 'distpar', [761.3, 0.2], 'dist', 'uniform');
settings.parameters(3) = struct('name', 'a111', 'distpar', [61.46, 0.2], 'dist', 'uniform');
settings.sample_points = 1e3;

output = saltelli_sobol(settings);