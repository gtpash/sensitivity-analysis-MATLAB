% test morris screening
%
% Author(s): Paul R. Miles | August 8, 2018

clear; close all; clc;

addpath('../')

c1 = 2;
c2 = 1;

settings.model = @(x) c1*x(1) + c2*x(2);

settings.parameters(1) = struct('name', 'Q1', 'mu', 0, 'sigma', 1, 'limits', [-1, 1]);
settings.parameters(2) = struct('name', 'Q2', 'mu', 0, 'sigma', 9, 'limits', [-2, 2]);
settings.sample_points = 200;
settings.sample_type = 'normal';
settings.delta = 1e-6;

% settings.parameters(1) = struct('name', 'Q1', 'mu', 0.5);
% settings.parameters(2) = struct('name', 'Q2', 'mu', 0.5);
% settings.sample_points = 200;
% settings.sample_type = 'uniform';

% rng(5)
output = morris_screening(settings)

% rmpath('../');