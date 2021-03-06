function parameters = setup_parameters(sample_type, user_parameters)
% function parameters = setup_parameters(sample_type, parameters)
%
% Description: Check to make sure sufficient requirements have been defined
% for the parameter structure given the sample type being used (normal or
% uniform).
%
% Author(s) Paul R. Miles | August 9, 2018

p = length(user_parameters);
mu = zeros(p,1);
for ii = 1:p
    mu(ii) = user_parameters(ii).mu;
end

defaults = default_parameters(mu, sample_type);

dp = fields(defaults);

for ii = 1:p
    for jj = 1:length(dp)
        if isfield(user_parameters, dp{jj})
            parameters(ii).(dp{jj}) = user_parameters(ii).(dp{jj});
        else
            parameters(ii).(dp{jj}) = defaults(ii).(dp{jj});
        end
    end
end