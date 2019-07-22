function B = generate_B(p)
% function b = generate_b(p)
%
% Generate the b matrix based on pg. 334 in [1].  The b matrix is the
% [(p + 1) x p] strictly lower triangular matrix of ones.
A = ones(p + 1, p);
B = tril(A, -1); % takes the lower triangular portion below the main diagonal