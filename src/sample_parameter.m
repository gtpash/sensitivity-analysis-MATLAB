function [samples] = sample_parameter(par,rsize,scale_grad)

% if not provided, default to false
if nargin <3
    scale_grad = false;
end

if scale_grad
    if strcmpi(par.dist, 'uniform')
        samples = rand(rsize);
    elseif strcmpi(par.dist, 'normal')
        samples = randn(rsize);
    else
        error('WARNING: Distribution type %s cannot be normalized\n',par.dist)
    end
else
    if strcmpi(par.dist,'normal')
        samples = normrnd(par.distpar(1),par.distpar(2),rsize);
    elseif strcmpi(par.dist,'uniform')
        a = par.distpar(1) - par.distpar(2)*abs(par.distpar(1));
        b = par.distpar(1) + par.distpar(2)*abs(par.distpar(1));
        samples = unifrnd(a,b,rsize);
    elseif strcmpi(par.dist,'beta')
        samples = betarnd(par.distpar(1),par.distpar(2),rsize);
    elseif strcmpi(par.dist,'gamma')
        samples = gamrnd(par.distpar(1),par.distpar(2),rsize);
    end
end