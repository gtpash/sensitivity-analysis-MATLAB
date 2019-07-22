% Test active subspace analysis

% setup workspace
clear; close all; clc;

%% Setup problem 
p = 2; % Number of parameters
num = 100; % Number of values to plot.
M = 100; % Gradient matrix sample number.
data.p = p;
data.M = M;
data.delta = 1e-6;
x1 = linspace(-1,1,num);
x2 = linspace(-1,1,num);

[X1,X2] = meshgrid(x1,x2);
x.x1 = X1;
x.x2 = X2;
% physical space sampling end points.
a = [-0.05;-0.05];
b = [0.05;0.05];
data.end_pts = [a b];

model = @(x) exp_fun(x);
Y = model(x);

%% Run Active Subspace Analysis
[asout] = active_subspace(model, data);

% %% Plot gradient map
% figure(1);
% surf(X1,X2,Y);
% colorbar
% 
% %% Plot response surface
% figure(2);
% plot(y_eval,f_approx_eval,'r');
% hold on
% plot(y_act,f_sam,'bo')