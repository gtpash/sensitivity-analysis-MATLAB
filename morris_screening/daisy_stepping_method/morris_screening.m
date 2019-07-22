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
% - levels: Integer, Grid levels
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
levels = settings.levels;
sample_points = settings.sample_points;

% determine number of parameters
p = length(parameters);

val = zeros(p,1);
for ii = 1:p
    val(ii) = parameters(ii).value;
end

Qstar = zeros(sample_points, p);
d = zeros(sample_points, p);
qstar = zeros(1, p);
for ii = 1:sample_points
    Num = rand(1, p);
    for jj = 1:levels - p
        for kk = 1:p
            v = abs(jj/levels - Num(kk));
            if v < val(kk)
                qstar(kk) = jj/levels;
                val(kk) = v;
            end
        end
    end
    Qstar(ii,:) = qstar;
    Delta = levels/(2*(levels-1));
    
    % generate matrices
    Bstar = generate_Bstar(qstar, Delta);
    
    % evaluate model
    y = zeros(p + 1, 1);
    for jj = 1:p + 1
        y(jj) = model(Bstar(jj, :));
    end
    
    % calculate elementary effect
    for jj = 1:p
        for kk = 1:p
            if Bstar(jj + 1, kk) ~= Bstar(jj, kk)
                d(ii,kk) = (y(jj+1) - y(jj))/Delta;
            end
        end
    end
end

% Construct the Morris indices.
mu = (1/sample_points)*sum(d, 1);
mu_i_star = (1/sample_points)*sum(abs(d), 1);
sigma_i = (1/(sample_points-1))*sum(bsxfun(@minus, d, mu).^2, 1);

% append features to output
output.mu_i_star = mu_i_star;
output.sigma_i = sigma_i;
output.Qstar = Qstar;
output.d = d;