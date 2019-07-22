function [dg] = dmapdtheta(parameters)
% function theta = dmapdtheta(theta0, parameters)
%
% Derivative of linear map that maps parameters from random sampling to the
% physical space for model evaluation.
%
% Author(s) Graham T. Pash | October 31, 2018
%           Paul R. Miles
%
%

% unpack parameters
p = length(parameters);
dg = zeros(p, 1);
for ii = 1:p
    if strcmpi(parameters(ii).dist, 'normal')
        dg(ii) = parameters(ii).distpar(2);
    elseif strcmpi(parameters(ii).dist, 'uniform')
        a = parameters(ii).distpar(1) - parameters(ii).distpar(2)*abs(parameters(ii).distpar(1));
        b = parameters(ii).distpar(1) + parameters(ii).distpar(2)*abs(parameters(ii).distpar(1));
        dg(ii) = (b - a);
    else
        error('WARNING: Unknown distribution type - %s\n', parameters(ii).type)
    end
end