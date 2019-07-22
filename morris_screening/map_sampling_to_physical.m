function theta = map_sampling_to_physical(theta0, parameters)
% function theta = map_sampling_to_physical(theta0, parameters)
%
% Map parameters from random sampling to the physical space for model
% evaluation.
%
% Author(s) Paul R. Miles | August 8, 2018
%
%

% unpack parameters
p = length(parameters);
theta = zeros(p, 1);
for ii = 1:p
    if strcmpi(parameters(ii).dist, 'normal')
        mu = parameters(ii).distpar(1);
        sigma = parameters(ii).distpar(2);
        theta(ii) = mu + sigma*theta0(ii);
    elseif strcmpi(parameters(ii).dist, 'uniform')
        a = parameters(ii).distpar(1) - parameters(ii).distpar(2)*abs(parameters(ii).distpar(1));
        b = parameters(ii).distpar(1) + parameters(ii).distpar(2)*abs(parameters(ii).distpar(1));
        theta(ii) = a + theta0(ii)*(b - a);
    else
        error('WARNING: Unknown distribution type - %s\n', parameters(ii).type)
    end
end