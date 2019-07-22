% test morris screening
%
% Author(s): Paul R. Miles | August 8, 2018

clear; close all; clc;

addpath('../')

P = linspace(0, 1);

settings.model = @(x) sum(x(1).*P.^2 + x(2).*P.^4 + x(3).*P.^6)/length(P);

settings.parameters(1) = struct('name', 'a1', 'value', 1);
settings.parameters(2) = struct('name', 'a11', 'value', 1);
settings.parameters(3) = struct('name', 'a111', 'value', 1);
settings.levels = 100;
settings.sample_points = 200;

% rng(5)
output = morris_screening(settings)