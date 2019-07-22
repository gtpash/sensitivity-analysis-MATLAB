function theta = map_sampling_to_physical(theta0, sample_type, parameters)
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
theta = zeros(1, p);
for ii = 1:p
    if strcmpi(sample_type, 'normal')
        theta(ii) = parameters(ii).mu + parameters(ii).sigma*theta0(ii);
    elseif strcmpi(sample_type, 'uniform')
        theta(ii) = parameters(ii).limits(1) + theta0(ii)*(parameters(ii).limits(2) - parameters(ii).limits(1));
    else
        error('WARNING: Unknown distribution type - %s\n', parameters(ii).type)
    end
end