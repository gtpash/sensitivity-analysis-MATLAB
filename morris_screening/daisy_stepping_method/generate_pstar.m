function Pstar = generate_pstar(p)
% function Pstar = generate_pstar(npar)
%
% Generate the P* matrix based on pg. 334 in [1].  The P* matrix is the
% [p x p] matrix constructed by randomly permuting the columns of a
% [p x p] identity matrix.

Pstar = eye(p); % generates [p x p] identity matrix
if rand(1) > 0.5 % permutes on average half the time
    Pstar = Pstar(:, end:-1:1);
end