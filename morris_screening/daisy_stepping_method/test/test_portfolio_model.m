% test morris screening
%
% Author(s): Paul R. Miles | August 8, 2018

clear; close all; clc;

addpath('../')

c1 = 2;
c2 = 1;

settings.model = @(x) c1*x(1) + c2*x(2);

settings.parameters(1) = struct('name', 'Q1', 'value', 1);
settings.parameters(2) = struct('name', 'Q2', 'value', 1);
settings.levels = 100;
settings.sample_points = 200;

% rng(5)
output = morris_screening(settings)