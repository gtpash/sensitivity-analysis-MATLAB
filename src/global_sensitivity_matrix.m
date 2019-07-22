function [output] = global_sensitivity_matrix(settings)
% function [output] = global_sensitivity_matrix(settings)
%
% Description: Run global sensitivity matrix to determine the most sensitive parameters
% in a model defined by the user.
%
% Author(s) Paul R. Miles | November 14, 2018
%           Graham T. Pash
%           Lider S. Leon
%           Alexandros Solomou
%
% Input: Structure
% - model: Anonymous model function
%       NOTE: model must accept a COLUMN vector as input
% - parameters: Cell array with names, nominal values
% - sample_points: Integer, number of sample points
%
% Output: Structure
% - 
%
% References:
% [1] Smith, R. C. (2013). Uncertainty quantification: theory, implementation,
%     and applications (Vol. 12). SIAM.
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

S = cell(sample_points,1);
n = length(h(samples(1,:)));
S_net = zeros(n, p);
mod_threshold = round(sample_points/10);
for jj = 1:sample_points
    if mod(jj,mod_threshold) == 0
        fprintf('Status: %i of %i\n', jj, sample_points);
    end
    
    % Construct D
    D = datasample([-delta,delta], p);

    % calculate elementary effect
    d = zeros(n,1);
    for ii = 1:p
        if strcmpi(use_method, 'cs')
            y = h(samples(jj,:)'+1i*D(ii)*I(:,ii));
            d = imag(y)/D(ii);
        else           
            y0 = h(samples(jj,:)');
            y = h(samples(jj,:)'+D(ii)*I(:,ii));
            d = (y - y0)/D(ii);
        end
        S{jj}(:,ii) = d(:);
    end
    S_net = S_net + S{jj};
end

% append features to output
output.d = d;
output.S = S;
output.Sbar = S_net./sample_points;
output.parameters = parameters;
output.samples=samples;