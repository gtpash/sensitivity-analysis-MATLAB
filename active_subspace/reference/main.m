clear all; clc; clf;
addpath('./PolyfitnTools/');

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

Y = exp_fun(x);

model_fun = @(x)exp_fun(x);

% Check that these 2 methods equivalent.
[U,S,G] = grad_mat(model_fun,data);
C = G*G';
[W,D] = eig(C);
D = rot90(D,2);
W = fliplr(W);
W_1 = -W(:,1); % the negative is for visualization - could switch back.
W_2 = W(:,2);
% Active variable and response surface
n_pts = 1e+2;
x_sam = (-1 + 2*rand(n_pts,p));
f_sam = model_fun(x_sam);
y_act = W_1'*x_sam';
poly_order = 3;
f_approx = polyfitn(y_act',f_sam,poly_order);
y_eval = linspace(min(y_act),max(y_act),num);
f_approx_eval = polyvaln(f_approx,y_eval);

% Plot gradient map
figure(1);
surf(X1,X2,Y);
colorbar

% Plot response surface
figure(2);
plot(y_eval,f_approx_eval,'r');
hold on
plot(y_act,f_sam,'bo')
