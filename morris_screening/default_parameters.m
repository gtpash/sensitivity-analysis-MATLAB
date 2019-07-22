function defaults = default_parameters(distpar, sample_type)

p = size(distpar,1);
delta = 0.20;

for ii = 1:p
    if strcmpi(sample_type{ii}, 'normal')
        defaults(ii) = struct('distpar', [0, 1], 'dist', 'normal');
    elseif strcmpi(sample_type{ii}, 'uniform')
        defaults(ii) = struct('distpar', [0, delta], 'dist', 'uniform');
    end
end