function [out] = active_subspace(model, data)
% Identify active parameter subspace for a given model
%
% Author(s):
%   Lider S. Leon
%   Paul R. Miles
%   Graham T. Pash
%
% Description:
%
% I/O:
%

%% Singular Value Decomposition
% Check that these 2 methods equivalent.
[U,S,G] = grad_mat(model, data);
C = G*G'; % C is symmetric and semi-positive definite by construction
[W,D] = eig(C);
D = rot90(D,2); % rotate array 90 degrees counterclockwise
W = fliplr(W); % flip order of columns in array

%% Define output structure
out.W = W;
out.U = U;
out.S = S;
out.G = G;
out.D = D;
% W_1 = -W(:,1);
% W_2 = W(:,2);

% %% Active variable and response surface
% n_pts = 1e+2; % number of sample points
% % Sample the training input values x_i with respect ot its probability
% % density function and construct corresponding responses q_i=h(x_i).
% xi = (-1 + 2*rand(n_pts,p));
% qi = model(xi);
% % Project the sampled values x_i onto the active subspace by using the
% % transformation y = W_1'*x_i
% y = W_1'*xi';
% z = W_2'*xi';
end

function [U,S,G] = grad_mat(model,data)

M = data.M;
p = data.p;
a = data.end_pts(:,1);
b = data.end_pts(:,2);
delta = data.delta;
D = datasample([-delta,delta],p);
T = diag(b-a);

G = zeros(p,M);

x_0 = rand(M,p);
theta_0 = (T*x_0' + a)';
y = model(theta_0);
I = eye(2);
for i = 1:M
    D = datasample([-delta,delta],p);
    for j = 1:p
        %%%%%%%%%%%
        %         % Use only if disproportional magnitudes present.
        %         x_diff = x_0(i,:) + D(j)*I(j,:);
        %         theta_diff = (T*x_diff' + a)';
        %%%%%%%%%%%
        theta_diff = theta_0(i,:) + D(j)*I(j,:);
        y_diff = model(theta_diff);
        G(j,i) = (y_diff - y(i))/D(j);
    end
end
[U,S] = svd(G/sqrt(M));
end