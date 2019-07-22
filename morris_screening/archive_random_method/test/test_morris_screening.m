% test morris screening
%
% Author(s): Paul R. Miles | August 8, 2018

clear; close all; clc;

addpath('../')

settings.model = @(x) test_function(x);
settings.parameters(1) = struct('name', 'val1', 'value', 1);
settings.parameters(2) = struct('name', 'val2', 'value', 1);
settings.levels = 100;
settings.sample_points = 200;

rng(5)
output = morris_screening(settings)