function Bstar = generate_Bstar(qstar, Delta)
% function Bstar = generate_Bstar(qstar)
%
% Generate the B* matrix based on pg. 334 in [1].  The B* matrix is
% comprised of p +  model realizations with the elements in the i-th row
% representing the parameter values used in the i-th evaluation.  This
% matrix is constructed so that for every column, j = 1, ..., p, there are
% two rows that differ only in their j-th component.  BY subtracting the
% elements inconsecutive rows, one can evaluate the p elementary effects
% associated with the random initial point q*.

p = length(qstar); % number of parameters

Pstar = generate_pstar(p);
Dstar = generate_dstar(p);
B = generate_B(p);
Jp1p = ones(p + 1, p);
Jp11 = ones(p + 1, 1);

Bstar = (Jp11*qstar + (Delta/2)*((2*B - Jp1p)*Dstar + Jp1p))*Pstar; % see pg. 334 in [1]

% check size of Bstar
[m, n] = size(Bstar);
if m ~= p + 1 || n ~= p
    error('Bstar should be [%i, %i] matrix.  Found Bstar as [%i, %i]', p + 1, p, m, n);
end