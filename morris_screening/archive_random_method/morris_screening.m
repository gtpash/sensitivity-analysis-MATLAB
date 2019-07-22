function [output] = morris_screening(settings)
% function [output] = morris_screening(settings)
%
% Description: Run morris screening to determine the most sensitive parameters
% in a model defined by the user.
%
% Author(s) Paul R. Miles | August 8, 2018
%
% Input: Structure
% - model: Anonymous model function
% - parameters: Cell array with names, nominal values
% - sample_points: Integer, number of sample points
%
% Output: Structure
% - mu_i_star
% - sigma_i
%
% References:
% [1] Smith, R. C. (2013). Uncertainty quantification: theory, implementation,
%     and applications (Vol. 12). SIAM.
% [2] Morris, M. D. (1991). Factorial sampling plans for preliminary
%     computational experiments. Technometrics, 33(2), 161-174.

% unpack input settings
model = settings.model;
parameters = settings.parameters;
sample_points = settings.sample_points;
sample_type = settings.sample_type;
delta = settings.delta;

parameters = setup_parameters(sample_type, parameters);

% determine number of parameters
p = length(parameters);
I = eye(p); % [p x p] identity matrix
d = zeros(p, 1);
G = zeros(sample_points, p);
for jj = 1:sample_points
    if mod(jj,5) == 0
        fprintf('Status: %i of %i\n', jj, sample_points);
    end
    % Construct D
    D = datasample([-1,1], p) * delta;
    % Generate set
    x0 = randn(p, 1);
    % map x to physical coordinates
    theta0 = map_sampling_to_physical(x0, sample_type, parameters);
    % evaluate model at physical coordinates
    y0 = model(theta0);
    % calculate elementary effect
    for ii = 1:p
       theta = map_sampling_to_physical(x0 + D(ii)*I(:,ii), sample_type, parameters);
       y = model(theta);
       d(ii) = (y - y0)/D(ii);
    end
    
    G(jj,:) = d(:)';
end

% Construct the Morris indices.
mu = (1/sample_points)*sum(G, 1);
mu_i_star = (1/sample_points)*sum(abs(G), 1);
sigma_i = (1/(sample_points-1))*sum(bsxfun(@minus, G, mu).^2, 1);

% append features to output
output.mu_i_star = mu_i_star;
output.sigma_i = sigma_i;
output.d = d;
output.G = G;
output.parameters = parameters;