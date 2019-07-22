function Dstar = generate_dstar(p)
% function Dstar = generate_dstar(npar)
%
% Generate the D* matrix based on pg. 334 in [1].  The D* matrix is the
% [p x p] diagonal matrix whose elements are randomly chosen from the set
% {-1, 1}.

Dstar = eye(p);
for jj=1:p
    if rand(1) < 0.5
        Dstar(jj,jj) = -1;
    else
        Dstar(jj,jj) = 1;
    end
end