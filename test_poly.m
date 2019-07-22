% test Example 15.17 from [1]
%
% Full test both estimating Sobol' indices and Morris sensitivity measures
% 
% Author(s): Graham T. Pash | November 1, 2018
%
% References
% [1] Smith, R. C. (2013). Uncertainty quantification: theory, implementation,
%     and applications (Vol. 12). SIAM.

addpath ./morris_screening/ ./saltelli_sobol/

settings.model = @prodmodel;

% All 6 Qi ~ U(0,1)
settings.parameters(1:6) = struct('distpar', [.5 1],'dist','uniform');
settings.sample_points = 1e3;

sobol = saltelli_sobol(settings);
morris = morris_screening(settings);

%% auxillary functions
function [Y] = prodmodel(x)
p = length(x);
a = [78 12 0.5 2 97 33];
Y = 1; % initialize product
for ii = 1:p
   gi = (abs(4*x(ii)-2)+a(ii))/(1+a(ii));
   Y = Y*gi;
end
end