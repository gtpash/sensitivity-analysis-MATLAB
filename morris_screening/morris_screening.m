function [output] = morris_screening(settings)
% function [output] = morris_screening(settings)
%
% Description: Run morris screening to determine the most sensitive parameters
% in a model defined by the user.
%
% Author(s) Paul R. Miles | August 8, 2018
%           Graham T. Pash
%
% updated: October 31, 2018
%
% Input: Structure
% - model: Anonymous model function
%       NOTE: model must accept a COLUMN vector as input
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
% [3] L Leon, K Coleman, N Bravo, R.C. Smith. Gradient Evaluations for
%     Models with Disproportional Parameter Values. August 9, 2018.

% unpack input settings
f = settings.model;
parameters = settings.parameters;
sample_points = settings.sample_points;

% error handling for optional parameters
if isfield(settings,'delta')
    delta = settings.delta;
else
    delta = 1e-6;
end
if isfield(settings,'observation_error')
    observation_error = settings.observation_error;
else
    observation_error = false;
end
if isfield(settings,'scale_grad')
    scale_grad = true;
else
    scale_grad = false;
end
if ~isfield(settings,'use_method')
    use_method = 'fd';
else
    use_method = settings.use_method;
end

parameters = setup_parameters(parameters);

% determine number of parameters
p = length(parameters);
I = eye(p); % [p x p] identity matrix
d = zeros(p, 1);
Ghat = zeros(sample_points, p);

% set up mapping functions
if scale_grad
    % Define scaling gradient functions
    g = @(x) map_sampling_to_physical(x, parameters);
    h = @(x) f(g(x));
else
    h = @(x) f(x);
end

% generate samples
samples = zeros(sample_points,p);
for pp = 1:p
    samples(:,pp) = sample_parameter(parameters(pp),[sample_points,1],scale_grad);
end

mod_threshold = round(sample_points/10);
for jj = 1:sample_points
    if mod(jj,mod_threshold) == 0
        fprintf('Status: %i of %i\n', jj, sample_points);
    end
    
    % Construct D
    D = datasample([-delta,delta], p);

    % only compute y0 once if using finite differences
    if ~strcmpi(use_method, 'cs')
        y0 = h(samples(jj,:)');
    end
    
    % calculate elementary effect
    for ii = 1:p
        if strcmpi(use_method, 'cs')
            y = h(samples(jj,:)'+1i*D(ii)*I(:,ii));
            d(ii) = imag(y)/D(ii);
        else
            y = h(samples(jj,:)'+D(ii)*I(:,ii));
            d(ii) = (y - y0)/D(ii);
        end
        
        % sigma-scale the elementary effect
        if observation_error ~= false
            if strcmpi(parameters(ii).dist, 'normal')
                d(ii) = d(ii)*sqrt(parameters(ii).distpar(2)/observation_error);
            end
        end
    end
    Ghat(jj,:) = d(:)';
end

% Construct the Morris indices.
mu = (1/sample_points)*sum(Ghat, 1);
mu_i_star = (1/sample_points)*sum(abs(Ghat), 1);
sigma_i = (1/(sample_points-1))*sum(bsxfun(@minus, Ghat, mu).^2, 1);

% append features to output
output.mu_i_star = mu_i_star;
output.sigma_i = sigma_i;
output.d = d;
output.G = Ghat;
output.parameters = parameters;