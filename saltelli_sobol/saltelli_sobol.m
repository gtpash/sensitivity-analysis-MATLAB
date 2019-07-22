function [out] = saltelli_sobol(settings)
% function [out] = saltelli_sobol(settings)
% 
% main script for generating Sobol indices using the Saltelli algorithm
% 
% Description: Run morris screening to determine the most sensitive parameters
% in a model defined by the user.
%
% Author(s) Graham T. Pash | October 24, 2018
%           Paul R. Miles
%
% updated: November 1, 2018
% 
% Input: Structure
% - model: Anonymous model function
% - parameters: Cell array with names, nominal values
% - sample_points: Integer, number of sample points
%
% Output: Structure
% - Si  : first order Sobol' indices
% - Ti : second order Sobol' indices
% 
% References
% [1] Smith, R. C. (2013). Uncertainty quantification: theory, implementation,
%     and applications (Vol. 12). SIAM.
% 
% [2] A. Saltelli, "Making best use of model evaluations to compute
%     sensitivity indices," Computer Physics Communications, 145, pp. 280-297, 2002.

% unpack settings
model = settings.model;
parameters = settings.parameters;
M = settings.sample_points;

p = length(parameters);
% preallocate space for parameter matrices and response vectors
A = zeros(M,p); B = zeros(M,p); C = zeros(M,p,p);
yA = zeros(M,1); yB = zeros(M,1); yC = zeros(M,p);

% construct A, B
for pp = 1:p
    A(:,pp) = sample_parameter(parameters(pp),[M,1]);
    B(:,pp) = sample_parameter(parameters(pp),[M,1]);
end

% construct Ci matrices
for pp = 1:p
    C(:,:,pp) = B;
    C(:,pp,pp) = A(:,pp);
end

mod_threshold = round(M/10);
% compute the model responses
for jj = 1:M
    if mod(jj,mod_threshold) == 0
        fprintf('Status: %i of %i\n', jj, M);
    end
    yA(jj) = model(A(jj,:));
    yB(jj) = model(B(jj,:));
    for pp = 1:p
        yC(jj,pp) = model(C(jj,:,pp));
    end
end

f02 = mean(yA)*mean(yB);
Si = zeros(p,1);
Ti = zeros(p,1);

% compute sensitivity indices
for pp = 1:p  
    Si(pp) = (((1/M)*yA'*yC(:,pp)) - f02)/((1/M)*(yA'*yA) - f02);
    Ti(pp) = 1-(((1/M)*yB'*yC(:,pp)) - f02)/((1/M)*(yA'*yA) - f02);
end

out.Si = Si;
out.Ti = Ti;