function [U,S,G] = grad_mat(model_fun,data)

M = data.M; % # of columns in gradient matrix
p = data.p; % # of parameters
a = data.end_pts(:,1); % lower parameter limits
b = data.end_pts(:,2); % upper parameter limits
delta = data.delta; % perturbation
% D = datasample([-delta,delta],p);
T = diag(b-a); % transformation matrix (mapping)

G = zeros(p,M); % allocate space

x_0 = rand(M,p); % 
theta_0 = (T*x_0' + a)';
y = model_fun(theta_0); % evaluate model at theta_0
I = eye(p);
for i = 1:M
    D = datasample([-delta,delta],p);    
    for j = 1:p
        %%%%%%%%%%%
%         % Use only if disproportional magnitudes present.
%         x_diff = x_0(i,:) + D(j)*I(j,:);
%         theta_diff = (T*x_diff' + a)';
        %%%%%%%%%%%
        theta_diff = theta_0(i,:) + D(j)*I(j,:);
        y_diff = model_fun(theta_diff);
        G(j,i) = (y_diff - y(i))/D(j);
    end
end
[U,S] = svd(G/sqrt(M));

