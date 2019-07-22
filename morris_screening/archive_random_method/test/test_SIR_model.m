% test SIR model


clear; close all; clc;

addpath('../')
addpath('SIR_codes/');



gamma = 0.2;
r = 0.6;
delta = 0.15;
k = 0.1;
params = [gamma r delta k];

S0 = 900;
I0 = 100;
R0 = 0;
N = 1000;

tf = 5;
dt = 0.1;
t_vals = 0:dt:tf;
Y0 = [S0; I0];
ode_options = odeset('RelTol',1e-6);

settings.model = @(x) SIR_model(N, t_vals, Y0, ode_options, x);

% [gamma r delta k];
settings.parameters(1) = struct('name', '\gamma', 'mu', gamma, 'limits', [0, 1]);
settings.parameters(2) = struct('name', 'k', 'mu', k, 'limits', [0, 1]);
settings.parameters(3) = struct('name', 'r', 'mu', r, 'limits', [0, 1]);
settings.parameters(4) = struct('name', '\delta', 'mu', delta, 'limits', [0, 1]);

settings.sample_points = 100;
settings.sample_type = 'uniform';
settings.delta = 1e-6;

% rng(5)
output = morris_screening(settings)