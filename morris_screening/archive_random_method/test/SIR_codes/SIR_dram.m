%   
%                  SIR_dram.m
%

  clear all
  close all

  global N Y0

%
%  Load data which consists of 61 points of infectious data.
%

%   load SIR_data

%
% Set parameters and initial conditions
%

  S0 = 900;
  I0 = 100;
  R0 = 0;
  N = 1000;

  gamma = 0.2;
  r = 0.6;
  delta = 0.15;
  k = 0.1;
  params = [gamma r delta k];

  tf = 6;
  dt = 0.1;
  t_vals = 0:dt:tf;
  Y0 = [S0; I0];

%
% Define the standard deviation sigma for the measurement error.
%

  sigma = 0.1*500;

%
% Integrate the system using ode45.m.  We exploit the fact that R(t) = N - S(t) - I(t).
%

  ode_options = odeset('RelTol',1e-6);
  [t,Y] = ode45(@SIR_rhs,t_vals,Y0,ode_options,params);
  Y(:,3) = N - Y(:,1) - Y(:,2);

  %error = sigma*randn(size(t_vals));
  %obs = Y(:,2) + error';

%
% Plot the results.
%

  figure(1)
  plot(t,Y(:,2),'x','linewidth',3)
  set(gca,'Fontsize',[22]);
  xlabel('Time')
  ylabel('Number of Infectious')