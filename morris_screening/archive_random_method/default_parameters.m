function defaults = default_parameters(mu, sample_type)

p = length(mu);
delta = 0.20;

if strcmpi(sample_type, 'normal')
    for ii = 1:p
        defaults(ii) = struct('mu', mu(ii), 'sigma', delta*abs(mu(ii)));
    end
else
    for ii = 1:p
        if mu(ii) == 0
            defaults(ii) = struct('mu', mu(ii), 'limits', [-delta, delta]);
        else
            defaults(ii) = struct('mu', mu(ii), 'limits', [mu(ii) - delta*abs(mu(ii)), mu(ii) + delta * abs(mu(ii))]);
        end
    end
end